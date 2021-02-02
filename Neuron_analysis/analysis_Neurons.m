%% Plot waveforms
% First load group TimeStamps and neuron TS
A = TS/10000;
n = 2; % Channel
C = zeros(length(TimeStamps),1);
for ii = 1:length(A)
    
    B = find(TimeStamps==A(ii,1));
    C(B) = 1;    
end

% Plot mean waveform
D = WaveForms(:,n,:);
E = D(logical(C),:,:);
F = mean(E(:,:));
plot(F);

% Plot all waveforms
G = E(:,:);
for ii = 1:length(E)
plot(G(ii,:),'k');
hold on;
end
hold on;
plot(F,'r');
hold off;

%% Waveform calculations
% First load group TimeStamps (GR1) and neuron TS (GR1_1)
% NOTE: Waveform time length = 900us (27units), 15 units = 500ms.
A = TS/10000;
n = 1; % Channel
C = zeros(length(TimeStamps),1);
for ii = 1:length(A)
    
    B = find(TimeStamps==A(ii,1));
    C(B) = 1;    
end

D = WaveForms(:,n,:); % Select the wave from n channel
E = D(logical(C),:,:);
allUP = E(:,:); 
allDOWN = -allUP; % Inverted values
mUP = mean(allUP);
mDOWN = -mUP; % Inverted values

% Only mean waveform
[pk2(1,1),pk2(2,1),w2,~] = findpeaks(mUP);
[pk,loc,~,pr] = findpeaks(mDOWN);
F = mDOWN>0;
I(1,1) = find(F,1,'first');
I(1,2) = find(F,1,'last');
    % Params
    wf.a = pr(1,1);
    wf.b = pr(1,2);
    wf.pk1(1,1) = pk(1,1);
    wf.pk1(1,2) = loc(1,1);
    wf.pk2 = transpose(pk2);
    wf.pk3(1,1) = pk(1,2);
    wf.pk3(1,2) = loc(1,2);
    wf.c = pk2(1,2) - pk1(1,2);
    wf.d = pk3(1,2) - pk2(1,2);
    wf.e = pk2(1,1) + pk1(1,1);
    wf.f = pk3(1,1) + pk2(1,1);
    wf.g = w2;
    wf.h = pk2(1,2) - I(1,1);
    wf.i = I(1,2) - pk2(1,2);

%% neuron ISI histogram
    
% Load TS from the neuron

% params
epsilon = 1e-100;
nBins = 500;

ISI = diff((TS)/10) + epsilon;

maxLogISI = max(real(log10(ISI)))+1;
minLogISI = floor(min(real(log10(ISI))));

H = ndhist(log10(ISI)', nBins, minLogISI, maxLogISI);
binsUsed = logspace(minLogISI,maxLogISI,nBins);

% Plot 
plot(binsUsed, H);
set(gca, 'XScale', 'log', 'XLim', [10^minLogISI 10^maxLogISI]);
set(gca, 'YTick', max(H));
hold on
plot([1 1], get(gca, 'YLim'), 'r:')
hold off

%% Align to ttl/ behaviour

% Before starting load the ttl var ant TimeStamps
% ttl = []; % ttl - min(ts); Array n x 2, first column is the time and second the behaviour. 
Dt = 408; % Time difference between the start of the video and the recording in sec

% Select group and neuron.
gr = 2; % Group number
nr = 1; % Neuron number
neuron = ['GR', num2str(gr),'_',num2str(nr), '.mat'];
TT = load(neuron,'TS');
TT = TT.TS;
TTb = ((TT(:,1)/10000)-min(ts)+Dt)/60; % For the video
% TTb = ((TT(:,1)/10000)+Dt); % For the shock

% Analysis of firing locked to behaviour
% secs = floor(ttl(:,1))*60 + (ttl(:,1)-floor(ttl(:,1)))/60*100; % Time in seconds
secs = ttl(:,1); 
figure; histogram(TTb,ceil(length(TTb)/2));
hold on;
stairs(secs,ttl(:,2)*10,'r');
xlim([min(secs) max(secs)]);
ylim([0 max(ttl(:,2)*10)+10]);
hold off;

% % Analysis of psth
% fs = 30000;
% pre = 5 *fs;
% post = 5 *fs;
% bin = 10000;
% [psth,psth2,ts1,ts2] = opto_psth(TTb*fs,ttl*fs,fs,bin,'pre',pre,...
% 'post',post,'chart',2);
% pre = 1 *fs;
% post = 1 *fs;
% bin = 1000;
% [psth,psth2,ts1,ts2] = opto_psth(TTb*fs,ttl*fs,fs,bin,'pre',pre,...
% 'post',post,'chart',2);

%% Get data from psth and save to mat for each neuron

% Variables to modify 
GR = 8; % Group of tetrodes to analyse
Dt = 0; % Delay time
ttl = TTL(:,1); % TTL you want to use
TTL0s = struct; % Valiable to be saved
sec = 0.5; % Seconds to analyse
name = ['TTL',num2str(floor(sec)),'s']; % Name of the valiable to be saved
bin = 1000;

tmp = dir(['GR',num2str(GR),'_','*.mat']); % all .mat that belong to the group
files = {tmp.name}'; 

if length(files)>9
    A = length(files)-9;
    for kk = 1:A
        B = files{1+kk};
        [files{1+kk}] = [];
        files{end+1} = B;
    end
    
    for kk = 1:length(files)
        file(kk) = ~isempty(files{kk});
    end
    files = files(file==1);
end

% Constant values
fs = 30000;
pre = sec *fs;
post = sec *fs;
excel = zeros(length(files),4);
Gr_N = cell(length(files),1);

for ii = 1:length(files)
    
    FstSpk_means = zeros(1,4);
    AllSpk_means = [];
    AllSpk = [];
    FstSpk = [];

    nr = ii;
    neuron = files{ii};
    TT = load(neuron,'TS');
    TT = TT.TS;
    TTb = ((TT(:,1)/10000)); %-(min(ts)+Dt)/60);
    
    % Calculate spike time around ttl
    [psth1,psth2,ts1,ts2] = opto_psth(TTb*fs,ttl*fs,fs,bin,'pre',pre,...
                    'post',post,'chart',2);
   
    % All spikes before & after ttl            
    for jj = 1:length(ts1)
        A = cell2mat(ts1(jj,1));
        if isempty(A) == 1
            A = NaN;
        end
        n = max(size(AllSpk,1),size(A,1));
        AllSpk(end+1:n,:) = nan;
        A(end+1:n,1) = nan;
        AllSpk = [AllSpk, A];
    end 
    AllSpk(:,1) = [];
    AllSpk = AllSpk/fs;
    
    % Spikes before & 1st after ttl 
    for kk = 1:length(ts2)
        ts3 = transpose(ts2);
        A = cell2mat(ts3(kk,1));
        if isempty(A) == 1
            A = NaN;
        end
        n = max(size(FstSpk,1),size(A,1));
        FstSpk(end+1:n,:) = nan;
        A(end+1:n,1) = nan;
        FstSpk = [FstSpk, A];
    end 
    FstSpk(:,1) = [];
    FstSpk = FstSpk/fs;
    
    % Find means of 1st spike after ttl
    FstSpk_means(1,1) = nanmean(FstSpk(FstSpk>0));
    FstSpk_means(1,2) = nanmedian(FstSpk(FstSpk>0));
    FstSpk_means(1,3) = nanstd(FstSpk(FstSpk>0));
    FstSpk_means(1,4) = length(FstSpk(FstSpk>0));
        
    % Calculate the firing peaks (and times) for all-spikes
    [psth_spx, psth_t] = psth_hist(psth1, bin);
    [AllSpk_means(:,1), AllSpk_means(:,2)] = findpeaks(psth_spx,...
            psth_t/fs,'MinPeakHeight',0.8);
    if isempty(AllSpk_means) == 1
        AllSpk_means(1,1) = NaN;
        AllSpk_means(1,2) = NaN;
    end
    
    TTL0s = struct('FstSpk', FstSpk, 'AllSpk', AllSpk,...
        'FstSpk_means', FstSpk_means, 'AllSpk_means', AllSpk_means);
    save(neuron,genvarname(name),'-append'); 
    saveas(gcf, genvarname(['TT',num2str(GR),'_',num2str(ii),'_',...
        num2str(floor(sec)),'s']), 'jpg')
    saveas(gcf, genvarname(['TT',num2str(GR),'_',num2str(ii),'_',...
        num2str(floor(sec)),'s']), 'svg')
    Gr_N{ii,:} = ['GR',num2str(GR),'_',num2str(ii)];
    excel(ii,1:4) = FstSpk_means;
    clearvars -except GR Dt ttl files fs pre post bin ts TTL...
        name sound excel sec Gr_N
%     close gcf;
end
% T = table('Size',[length(files) 5],'VariableTypes',...
%     {'string','double','double','double','double'},...
%     'VariableNames',{'Neuron','Mean','Median','Std','numSpk'});
% T.Neuron = Gr_N;
% T.Mean = excel(:,1);
% T.Median = excel(:,2);
% T.Std = excel(:,3);
% T.numSpk = excel(:,4);
% 
% writetable(T, ['TT',num2str(GR),'shock.xls']);
