function [wdata, wfreq] = wavelet_analysis(signal,time,tsegment,freq,varargin)

% The continuous wavelet transform analysis is similar to the FFT, but
% shows the temporal evolution of the dominant frequency.
%
% INPUTS: 
%   - signal: signal to analyse, a vector
%   - time: time corresponding to the read signal, an interval: [ts1 ts2]
%   - tsegment: time period to analyse, an interval: [t1 t2]
%   - freq: frequency range to analyse, an interval: [f1 f2]
%   Varargins
%   - fr: sampling rate. Use as 'fr',value. Default 30000.
%   - dj: something like the resoution. Use as 'dj',value. Default 0.08.
%   - n:  downsampling factor. Use as 'n',value. Default 3. 
%   - mother: the mother wavelet function. Use as 'mother','name'. Default
%             'Morlet'. Also choices 'PAUL', or 'DOG'.
%   - tInt: time intervals to segment before analysing, in seconds. Use as
%           'tInt',value. NO default value, all time will be analysed.
%
% OUTPUTS: 
%   - wdata: the WAVELET transform of signal.
%   - wfreq: The frequency range analysed.
% Results are automatically saved to .dat files named:
%   - 'signal_wavelet.dat': is the WAVELET transform of signal. FLOAT(WAVE)
%     gives the WAVELET amplitude, ATAN(IMAGINARY(WAVE),FLOAT(WAVE) gives
%     the WAVELET phase. The WAVELET power spectrum is ABS(WAVE)^2. Its
%     units are sigma^2 (the time series variance).
%   - 'signal_freq.dat': The frequency range analysed.
% To read the saved .dat please refer to the bottom of this function.
%
% The coi is calculated. If you want to plot it just uncomment the
% corresponding section.
% 
% Examples: 
% To only save the data
% wavelet_analysis(signal,time,[100 300],[0.5 300]);
% wavelet_analysis(signal,time,[100 300],[0.5 300],'fr',42,'dj',42,'n',42,...
%                  'mother','Morlet','tInt',42);
% To keep the result in the workspace
% [wdata, wfreq] = wavelet_analysis(signal,time,[100 300],[0.5 300]);
%
% In case of "out of memory" error, try including a tInt and pre-segment
% the analysis.
% 
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2019
% Laboratory of Network Neurophysiology
% Institute of Experimantal Medicine, Hungary.
%
% Wavelet software was provided by C. Torrence and G. Compo. See wavelet
% function for more details.
%
% MATLAB toolboxes: Signal Processing Toolbox, Statistics and Machine
% Learning Toolbox.
% -------------------------------------------------------------------------

% Optional input variables ------------------------------------------------
% Params, defalut values:
fr = 30000;
dj = 0.08; % 0.05, As set by Hanya B.    
n = 3;
mother = 'Morlet';

if nargin
  for i=1:2:size(varargin,2)
    switch varargin{i}
      case 'fr'
        fr = varargin{i+1};
      case 'dj'
        dj = varargin{i+1};
      case 'n'
        n = varargin{i+1};
      case 'mother'
        mother = varargin{i+1};
      case 'tInt'
        tInt = varargin{i+1};
      otherwise
        errordlg('Unknown varargin');
    end
  end
else
end
%--------------------------------------------------------------------------

% Generate the files to be saved
out = inputname(1);
wave = [out,'_',num2str(tsegment(1)),'_',num2str(tsegment(2)),'_wavelet.dat'];
f = [out,'_',num2str(tsegment(1)),'_',num2str(tsegment(2)),'_freq.dat']; 

% Keep definig params
dt = 1/fr/n; % amount of time between each Y value, i.e. the sampling time.
p0 = 1/freq(1);
pn = 1/freq(2);
s0 = pn/1.03; % the smallest scale of the wavelet.  Default is 2*DT.
sn = p0/1.03;
pad = 1;
j1 = ceil(log2(sn/s0)/dj); %the # of scales minus one. Scales range from S0 up to S0*2^(J1*DJ),
%        to give a total of (J1+1) scales. Default is J1 = (LOG2(N DT/S0))/DJ.

%     % Wavelet params as set by Hanya B
%     dt = 1 / fr/n;
%     pad = 1;  
%     j1 = ceil((1/dj) * log2(length(signal)/2));
%     j1 = ceil(j1);
%     s0 = 2 * dt; 

t0 = tsegment(1);
tf = tsegment(2);
if exist('tInt','var')==1
    % Make the time segments
    tArray = [];
    t1 = t0; t2 = t1+tInt;
    tArray = [tArray; [t1 t2]];
    while t2 <= (tf-1)
    t1 = t2; t2 = t1+tInt; % With overlap t1 = t2-10; t2 = t1+tInt;
    tArray = [tArray; [t1 t2]];
    end
elseif exist('tInt','var')==0
    tArray = [t0 tf];
end

for ii = 1:size(tArray,1)
    if ii == 1
        int = [tArray(ii,1) tArray(ii,2)]; 
        in = (int(1)-time(1))*(fr/n)+1;
        fin = (int(2)-time(1))*(fr/n);
        S = downsample(signal(1,in:fin),n);
        
        [wt,period,~,coi] = wavelet(S,dt,pad,dj,s0,j1,mother);
        power = zscore((abs(wt)).^2);
        fq = 1 ./ period;
        % l = log(wpow) - log(repmat(mean(wpow,2),1,size(wpow,2))); % wpow
        % = matriz de la potencia del wavelet.

        figure;
        surface('XData',[int(1) int(2)],'YData',[max(fq) min(fq)],...
            'CData',power,'ZData', zeros(2,2), 'CDataMapping','scaled',...
            'FaceColor','texturemap', 'EdgeColor', 'none');
%         xlim([int(1) int(2)]);
        set(gca, 'YScale', 'log');
        ylim([min(fq) max(fq)]);
        colormap jet
        colorbar; ylabel('Frequency (Hz)'); 
%         hold on; % Uncomment to plot the coi.
%         a = area(linspace(int(1),int(2),length(coi)),log2(abs(coi))); 
%         a.FaceColor = [0.6 0.6 0.6];
%         a.FaceAlpha = 0.8;
%         a.EdgeColor = [0.6 0.6 0.6];
%         hold off;
        
        % Saving into .dat
        % initialize a new .dat file
        fid = fopen(wave,'w+');    
        wt_real = real(wt);
        wt_imag = imag(wt);
        wt_r1i = [wt_real wt_imag];
        fwrite(fid,wt_r1i,'int16'); % Remember the data type to open it.
        fclose(fid);
        
        fid = fopen(f,'w+');    
        fwrite(fid,fq(:,:),'int16'); % Remember the data type to open it.
        fclose(fid);
    
        clear wt coi power wt_imag wt_real

    elseif ii < length(tArray)
        int = [tArray(ii,1) tArray(ii,2)]; 
        in = (int(1)-time(1))*(fr/n)+1;
        fin = (int(2)-time(1))*(fr/n);
        S = downsample(signal(1,in:fin),n);
        
        [wt,period,~,~] = wavelet(S,dt,pad,dj,s0,j1,mother);
        power = zscore((abs(wt)).^2);
        fq = 1 ./ period;

        surface('XData',[int(1) int(2)],'YData',[max(fq) min(fq)],...
            'CData',power,'ZData', zeros(2,2), 'CDataMapping','scaled',...
            'FaceColor','texturemap', 'EdgeColor', 'none');
        colormap jet
        colorbar;
%         hold on; % Uncomment to plot the coi.
%         a = area(linspace(int(1),int(2),length(coi)),log2(abs(coi))); 
%         a.FaceColor = [0.6 0.6 0.6];
%         a.FaceAlpha = 0.8;
%         a.EdgeColor = [0.6 0.6 0.6];
%         hold off;
        
        % Saving into .dat
        fid = fopen(wave,'a+');
        wt_real = real(wt);
        wt_imag = imag(wt);
        wt_r1i = [wt_real wt_imag];
        fwrite(fid,wt_r1i,'int16'); % Remember the data type to open it.
        fclose(fid);
        
        fid = fopen(f,'a+');    
        fwrite(fid,fq(:,:),'int16'); % Remember the data type to open it.
        fclose(fid);  
       
        clear wt coi power wt_imag wt_real
        
    elseif ii == length(tArray)
        int = [tArray(ii,1) tf]; 
        in = (int(1)-time(1))*(fr/n)+1;
        fin = (int(2)-time(1))*(fr/n);
        S = downsample(signal(1,in:fin),n);
        
        [wt,period,~,coi] = wavelet(S,dt,pad,dj,s0,j1,mother);
        power = zscore((abs(wt)).^2); 
        fq = 1 ./ period;

        surface('XData',[int(1) int(2)],'YData',[max(fq) min(fq)],...
            'CData',power,'ZData', zeros(2,2), 'CDataMapping','scaled',...
            'FaceColor','texturemap', 'EdgeColor', 'none');
        colormap jet
        colorbar;
%         hold on; % Uncomment to plot coi
%         N = length(coi); % Last plot
%         t3 = linspace(int(1),int(2),N);
%         coi2 = transpose(coi(N/n+1:end,1));
%         a = area(t3,coi2);
%         a.FaceColor = [0.6 0.6 0.6];
%         a.FaceAlpha = 0.8;
%         a.EdgeColor = [0.6 0.6 0.6];
%         hold off;
        
        % Saving into .dat
        fid = fopen(wave,'a+');
        wt_real = real(wt);
        wt_imag = imag(wt);
        wt_r1i = [wt_real wt_imag];
        fwrite(fid,wt_r1i,'int16'); % Remember the data type to open it.
        fclose(fid);
        
        fid = fopen(f,'a+');    
        fwrite(fid,fq(:,:),'int16'); % Remember the data type to open it.
        fclose(fid);
        
        clear wt coi power wt_imag wt_real

    end
end

if nargout > 0
    wdata = wt_r1i;
    wfreq = fq(:,:);
else
end

end

% % To read the saved .dat, change the Data and time variables
% % Fist open the freq.dat to know number of rows;
% Data = 'mT1_300_600'; % Where 300 and 600 is the time
% time = [300 600];
% 
% file = dir(strcat(Data, '_freq.dat'));
% fid = fopen(file.name, 'r');
% freq = fread(fid, 'int16'); 
% fclose(fid);
% % freq is a variable of r rows and 1 col.
% % Then read the bytes of the wave
% file = dir(strcat(Data, '_wavelet.dat'));
% fid = fopen(file.name,'r');
% r = size(freq,1);
% same_real = fread(fid, [r file.bytes/r/2/2], 'int16'); % Where r is the number  
% % of rows as read in freq, 2 is the bytes number in int16 and 2 is because
% % the data is double (real and imaginary)
% same_imag = fread(fid, [r file.bytes/r/2/2], 'int16'); % Where r is the number  
% % of rows as read in freq, 2 is the bytes number in int16 and 2 is because
% % the data is double (real and imaginary)
% fclose(fid);
% wave = complex(same_real, same_imag);
% power = zscore((abs(wave)).^2);
% surface('XData',[time(1) time(2)],'YData',[max(freq) min(freq)],...
%     'CData',power,'ZData', zeros(2,2), 'CDataMapping','scaled',...
%     'FaceColor','texturemap', 'EdgeColor', 'none');
