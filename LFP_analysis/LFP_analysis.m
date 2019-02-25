%% Convert .continuous to dat 
% Important: In case of multiple .continuous for each channel. File name:
% XXX_CHn, where XXX is the number of processor and c the number of
% channel.
% Look into the settings.xml file. Locate the processor corresponding to
% the "source" (RhythmFPGA) and load only these files.
% Careful: Big numser of channel may bring 'out of memory'.

CH = dir('100_CH*.continuous'); % Loads 100 processor files.
nCH = length(CH);
fd = CH.folder;

% Depends on the geometry and type of electrode. For tetrodes:
nT = nCH/4;
n = 0;
for ii = 1:nT
    T.data = [];
    T.timestamps = [];
        
    for jj = 1:4
        ch = jj+n;
        chN = ['100_CH',num2str(ch),'.continuous'];
        [data, timestamps, info] = load_open_ephys_data([fd,'\',chN]);
        T.data(ch,:) = data;
        T.ts(1,:) = timestamps;
    end
    
    data = ['T',num2str(ii),'_data.dat'];
    %initialize a new .dat file
    fid = fopen(data,'w+');    
    % write data to raw .dat file
    fwrite(fid,T.data(:,:),'int16');
    % close dat file
    fclose(fid);
    
    ts = ['T',num2str(ii),'_ts.dat'];
    %initialize a new .dat file
    fid = fopen(ts,'w+');    
    % write data to raw .dat file
    fwrite(fid,T.ts(:,:),'int16');
    % close dat file
    fclose(fid);
    
    disp(['Files ',data,' and ',ts,' converted.']);
    clear T
    n = n+4;
end

%% Read .dat to MATLAB
% Implemented for tetrodes (num_channels = 4). Reads specific periods
% selected with the time variable.

TT = 1;
num_channels = 4;
fr = 30000;
time = [800 1400];

file = dir(['T',num2str(TT),'_data.dat']);
T = genvarname(['T',num2str(TT)]);
num_samples = file.bytes/(num_channels * 2); % uint16 = 2 bytes
Mt = num_samples/fr;
t = linspace(0,Mt,num_samples);
t2 = t(:,(time(1)*fr:time(2)*fr));

fid = fopen(file.name, 'r');
ftell(fid); % This sets the starting reading value to 0. The 'ans' should be 0.
fseek(fid,num_channels*time(1)*fr*2,'bof'); % Sets the starting reading at the begining of the part you are interested in (*2 for int16) 
Int = (time(2)*fr - time(1)*fr)+1; % Sets the interval to read.
ftell(fid);
rT = fread(fid, [num_channels, Int], 'int16'); %Reads the specified interval within the file.
fclose(fid);
eval ([T '= rT * 0.195;']);
clear rT

%% Wavelet using wavelet toolbox within MATLAB

% Params: Change these values to adjust them to your needs.
fr = 30000; % Sampling rate
int = [500 540];
in = ((int(1)-500)*fr)+1;
fin = (int(2)-500)*fr;
freq = [0.5 300]; % frequency limits
t = time2; % Time vector: < num_channels = 4; fileinfo = dir('T1.dat'); 
        % num_samples = fileinfo.bytes/(num_channels * 2); % uint16 = 2
        % bytes; Mt = num_samples/30000; t = linspace(0,Mt,num_samples);
signal = mT(1,in:fin); % The variable to analyse
mother = 'amor'; % The wavelet to use: 'amor' makes te Morlet
n = 3; % Downsampling factor
signal2 = downsample(signal,n);

[wt,f,coi] = cwt(signal2,mother,fr/n,'FrequencyLimits',freq);
power = zscore((abs(wt)).^2); 


% figure; subplot(2,1,1);
% surface('XData',[min(t) max(t)],'YData',[max(f) min(f)],'CData',abs(wt),...
%     'ZData', zeros(2,2), 'CDataMapping','scaled', 'FaceColor',...
%     'texturemap', 'EdgeColor', 'none');
% % xlim([min(t) max(t)]);
% xlim([int(1) int(2)]);
% set(gca, 'YScale', 'log');
% ylim([min(f) max(f)]);
% colormap jet
% colorbar; title('Wavelet spectrogram'); ylabel('Frequency (Hz)'); 
% subplot(2,1,2);
surface('XData',[min(t) max(t)],'YData',[max(f) min(f)],'CData',power,...
    'ZData', zeros(2,2), 'CDataMapping','scaled', 'FaceColor',...
    'texturemap', 'EdgeColor', 'none');
% xlim([min(t) max(t)]);
xlim([int(1) int(2)]);
set(gca, 'YScale', 'log');
ylim([min(f) max(f)]);
colormap jet
colorbar; ylabel('Frequency (Hz)'); xlabel('t (s)');
hold on; a = area(downsample(t(1,in:fin),n),transpose(coi)); % Only first plot
a.FaceColor = [0.6 0.6 0.6];
a.FaceAlpha= 0.8;
a.EdgeColor= [0.6 0.6 0.6];
hold off;


%% Wavelet using Torrence, C. and G. P. Compo, 1998. 
% Implemeted as new function: wavelet_analysis

fr = 30000; % Sampling rate
freq = [0.5 300];
dj = 0.1;
int = [700 800];
in = (int(1)-500)*fr;
fin = (int(2)-500)*fr;
signal = downsample(mT(1,in:fin),3);

dt = 1/fr/3;
p0 = 1/freq(1);
pn = 1/freq(2);
s0 = pn/1.03;
sn = p0/1.03;
pad = 1;
j1 = log2(sn/s0)/dj;
dj = dj;
mother = 'Morlet';
    
[wt,period,scale,coi] = wavelet(signal,dt,pad,dj,s0,j1,mother);
power = zscore((abs(wt)).^2); 

figure; 
% subplot(2,1,1);
% surface('XData',[min(t) max(t)],'YData',[max(freq) min(freq)],'CData',abs(wt),...
%     'ZData', zeros(2,2), 'CDataMapping','scaled', 'FaceColor',...
%     'texturemap', 'EdgeColor', 'none');
% xlim([min(t) max(t)]);
% set(gca, 'YScale', 'log');
% ylim([min(freq) max(freq)]);
% colormap jet
% colorbar; title('Wavelet spectrogram'); ylabel('Frequency (Hz)'); 
% subplot(2,1,2);
surface('XData',[min(t) max(t)],'YData',[max(freq) min(freq)],'CData',power,...
    'ZData', zeros(2,2), 'CDataMapping','scaled', 'FaceColor',...
    'texturemap', 'EdgeColor', 'none');
% xlim([min(t) max(t)]);
xlim([int(1) int(2)]);
set(gca, 'YScale', 'log');
ylim([min(freq) max(freq)]);
colormap jet
colorbar; ylabel('Frequency (Hz)'); xlabel('t (s)');

%% Fourier PS 
% Implemeted as new function: fourier_analysis

% with butterworth filter
% signal = signal2;
signal = T1(1,1:3000000);
fr = 30000; % Sampling rate
freq = [0.5 300];

wn = freq/(fr/2);   %bandpass
[b,a] = butter(3,wn);  
ft = filter(b,a,signal);
N = 16384; % 4096 % Higher values make more precise pics.
fq = linspace(0,fr,N);
F = fft(ft,N);
maxFreq = N/8; %~2756 Hz.
figure; 
plot(fq(1:maxFreq),abs(F(1:maxFreq)));
xlim([min(freq) max(freq)]);

% Multi-taper fourier transform - continuous data 
n = 3; % Downsampling factor
signal2 = downsample(signal,n);
params.tapers = [3 5];
params.Fs = fr/n;
params.fpass = freq;

[F,f]=mtspectrumc(signal2,params);
figure;
plot(f,smooth(F,1,'lowess')); % Smooth level 0.02 is good usually
xlim([min(freq) max(freq)]);

