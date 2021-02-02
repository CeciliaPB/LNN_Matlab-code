
% Params
Tab = Neurons_List;
fs = 30000;
psth_bins = 600; % 600 = 20ms
Stim = 'SK';

% Setting folders
mainFolder = pwd;
folder = Tab.Recording;

% Variable to save
PSTHall = zeros(size(Tab,2),ceil(36001/psth_bins));

for kk = 1:size(folder,1)
    
    % Go to folder
    cd([mainFolder filesep folder{kk,:}]);
    
    % Load TTLs & Choose variable to use as stimulus
    if kk == 1 
        load('TTLs.mat');
        switch Stim
            case 'SK'
                ttl = TTL; % The SK TTL is already called TTL. 
            case 'WN'
                if exist('NewWN','var')
                    ttl = NewWN;
                elseif exist('NewWS','var')
                    ttl = NewWS;
                elseif exist('TTL_WN','var')
                    ttl = TTL_WN;
                elseif exist('TTL_WN_video','var')
                    ttl = TTL_WN_video;
                elseif exist('sound','var')
                    ttl = sound;
                end
            case 'PT'    
                if exist('NewPT','var')
                    ttl = NewPT;
                elseif exist('TTL_PT','var')
                    ttl = TTL_PT;
                elseif exist('TTL_PT_video','var')
                    ttl = TTL_PT_video;
                elseif exist('sound','var')
                    ttl = sound;
                end
        end
    elseif min(folder{kk,:} == folder{kk-1,:}) == 1
    elseif min(folder{kk,:} == folder{kk-1,:}) == 0
        clearvars ttl TTL NewWN NewWS TTL_WN TTL_WN_video NewPT TTL_PT...
            TTL_PT_video sound
        load('TTLs.mat');
        switch Stim
            case 'SK'
                ttl = TTL; % The SK TTL is already called TTL. 
            case 'WN'
                if exist('NewWN','var')
                    ttl = NewWN;
                elseif exist('NewWS','var')
                    ttl = NewWS;
                elseif exist('TTL_WN','var')
                    ttl = TTL_WN;
                elseif exist('TTL_WN_video','var')
                    ttl = TTL_WN_video;
                elseif exist('sound','var')
                    ttl = sound;
                end
            case 'PT'    
                if exist('NewPT','var')
                    ttl = NewPT;
                elseif exist('TTL_PT','var')
                    ttl = TTL_PT;
                elseif exist('TTL_PT_video','var')
                    ttl = TTL_PT_video;
                elseif exist('sound','var')
                    ttl = sound;
                end
        end
    end
    
    % Find neuron
    NeuronID = dir(['*',num2str(Tab.Group(kk)),'_',num2str(Tab.Neuron(kk)),'.mat']);
    NeuronID = NeuronID.name;

    load(NeuronID,'TS')
    TT = ((TS(:,1)/10000));

    [psth1, ts1, psth1st, ts1st] = ttl_psth (TT*fs, ttl*fs, psth_bins,...
        'pre', 0.5, 'post', 0.7);
    [psth_spx, psth_t] = psth_hist(psth1, psth_bins);
%     hold on;
%     plot(psth_t/fs,psth_spx)
%     hold off;

    PSTHall(kk,:) = transpose(zscore(psth_spx));

end

    PSTHt = psth_t;

%%
ifr_bins = 100; % Sampling 0.1s
iFRn = [];
for ii = 1:length(TTL)
time = [TTL(ii)-0.5,TTL(ii)+0.5];
timestamps = TT(TT > time(1) & TT < time(2));
[iFR,~] = instFiringRate(timestamps,time,bins);
iFRn(ii,:) = iFR;
end
t = linspace(-0.5,0.5,size(iFRn,2));
figure;plot(t,iFRn)
iFRm = mean(iFRn);
hold on;
plot(t,iFRm,'k','Linewidth',2)

