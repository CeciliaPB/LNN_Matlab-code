function psth = allneurons_psth_old_version(ngroup, ttl, varargin)
% Get data from psth and save to mat for each neuron
% IMPORTANT: For the analysis of pairs use ttl_psth.
% 
% INPUTS: 
%   - ngroup: groups of neurons to perform the analysis. It can be a single
%   group (ngroup = 2, for GR2) or several (ngroup = [1:8], for GR1 to
%   GR8). 
%   - ttl: a vector containing the ttl timestamps.
% Varargin
%   - 'fs': sampling rate. Dafault 30000.
%   - 'Dt': Delay time, between the recording and the ttl. Ex. in case of
%   getting the ttl from a video. Default NOT used.
%   - 'pre': time before trigger to include in psth. Default 1s. 
%   - 'post': time after trigger to include in psth. Default 1s.
%   - 'bin': width of the histogram column. Adjust the 'bin', 1000 bin for
%   0.5s; 3000 bin for 2s. 
%   - 'excel': to save an excel of the psth. Default NOT save excel 
%   - 's': to save the images. Default NOT save figs
% 
% OUTPUTS: 
%   - psth: a structure with the spikes around the event (all and up to
%   1st), and the mean values for the firing times of the neuron.
%   
% Examples: 
% psth = allneurons_psth([1:8], ttl);
% psth = allneurons_psth(2, ttl, 'pre',0.5,'post',0.5,'fr',30000);
% 
% -------------------------------------------------------------------------
% Cec?lia Pardo-Bellver, 2020
% Laboratory of Network Neurophysiology
% Instinute of Experimantal Medicine, Hungary.
%
% MATLAB toolboxes: - 
% -------------------------------------------------------------------------
  
% Default params
fs = 30000;
pre = 0.05;
post = 0.05;
bin = 100;
excel = 0; % Default NOT save excel
s = 0; % Default NOT save figs
if nargin
    for ii = 1:length(varargin)
        switch varargin{ii}
            case 'fs'
                fs = varargin{ii+1};
            case 'Dt'
                Dt = varargin{ii+1};
            case 'pre'
                pre = varargin{ii+1};
                if pre >= 5000
                elseif pre < 5000
                    pre = pre * fs;
                end
            case 'post'
                post = varargin{ii+1};
                if post >= 5000
                elseif post < 5000
                    post = post * fs;
                end
            case 'bin'
                bin = varargin{ii+1};
            case 'excel'
                excel = 1;                
            case 's'
                s = 1;
                            otherwise
        end
    end
end
for nn = ngroup
tmp = dir(['GR',num2str(nn),'_','*.mat']); % all .mat that belong to the group
files = {tmp.name}'; 
% This reorders the files so after 1 comes 2 and not 10. Thank you MATLAB.
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
% Preallocate vars
psth = struct;
xls = zeros(length(files),4); % To save an excel with relevant values
Gr_N = cell(length(files),1);
name = ['psth',num2str(floor(pre/fs)),'s']; % Valiable to be saved
for ii = 1:length(files)
    
    FstSpk_means = zeros(1,4);
    AllSpk_means = [];
    AllSpk = [];
    FstSpk = [];
    neuron = files{ii};
    load(neuron,'TS');
      if exist('Dt', 'var') == 0
        TT = ((TS(:,1)/10000)); 
    elseif exist('Dt', 'var') == 1
        A = dir('*.continuous');
        B = char({A.name});
        [~, ts, ~] = load_open_ephys_data(B(1,:));
        TT = ((TS(:,1)/10000)-(min(ts)+Dt)/60);
    end
    
    % Calculate spike time around ttl
    [psth1,ts1,~,ts2] = ttl_psth(TT*fs, ttl*fs, bin, 'fs', fs, 'pre',...
        pre, 'post', post, 'chart', 2);
        sgtitle(['NeuronID: ' neuron(1:end-4)],'Interpreter','none');
   
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
        A = cell2mat(ts2(kk,1));
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
    [AllSpk_means(:,1), AllSpk_means(:,2)] = findpeaks(psth_spx,psth_t/fs,...
            'MinPeakHeight',0.8);
    if isempty(AllSpk_means) == 1
        AllSpk_means(1,1) = NaN;
        AllSpk_means(1,2) = NaN;
    end
    
    psth = struct('FstSpk', FstSpk, 'AllSpk', AllSpk,...
        'FstSpk_means', FstSpk_means, 'AllSpk_means', AllSpk_means);
    eval([name '= psth']);
    save(neuron,name,'-append');
      if s == 1
        saveas(gcf, matlab.lang.makeValidName([neuron(1:end-4),...
            '_', num2str(floor(post/fs)),'s']), 'jpg')
%         saveas(gcf, matlab.lang.makeValidName([group(1:end-4),'_',...
%             num2str(ii),'_',num2str(floor(post/fs)),'s']), 'svg')
    else
    end
    
    Gr_N{ii,:} = neuron(1:end-4);
    xls(ii,1:4) = FstSpk_means;
    clearvars -except GR Dt ttl files fs pre post bin ts TTL...
        name sound excel sec Gr_N s excel group psth
%     close gcf;
end
if excel == 1
    T = table('Size',[length(files) 5],'VariableTypes',...
        {'string','double','double','double','double'},...
        'VariableNames',{'Neuron','Mean','Median','Std','numSpk'});
    T.Neuron = Gr_N;
    T.Mean = xls(:,1);
    T.Median = xls(:,2);
    T.Std = xls(:,3);
    T.numSpk = xls(:,4);
    writetable(T, [neuron(1:end-6),'.xls']);
    
else
end
end
end