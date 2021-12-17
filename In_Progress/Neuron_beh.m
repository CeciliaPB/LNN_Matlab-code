
% see behaviourly driven cells.
% A lot of work to do here.
%
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2018
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
% -------------------------------------------------------------------------

for xx = 1:8
    
GR = xx; % Group of tetrodes to analyse
Dt = 1334;

TT = load(['GR',num2str(GR),'.mat']);
tmp = dir(['GR',num2str(GR),'_','*.mat']); % all .mat of the group
files = {tmp.name}'; 

for ii = 1:length(files)

    nr = ii;
    neuron = files{ii};
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
%     xlim([min(secs) max(secs)]);
    ylim([0 max(ttl(:,2)*10)+10]);
    hold off;
    
    saveas(gcf, genvarname(['GR',num2str(GR),'_',num2str(nr),...
        '_Beh']), 'svg');
    saveas(gcf, genvarname(['GR',num2str(GR),'_',num2str(nr),...
        '_Beh']), 'jpg');
end
end