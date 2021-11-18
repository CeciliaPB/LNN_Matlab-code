function [TTL, pulseon, pulseoff] = load_events(sessionpath)

%   Creates a structure of events, each eventChannel in one variable.
%   INPUTS:
%       - sessionpath: string, folder where the all_channels.events is
%                      located.       
% -------------------------------------------------------------------------
% Cec�lia Pardo-Bellver, 2021
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
% -------------------------------------------------------------------------

% Load open ephys events file
[data, timestamps, info] = load_open_ephys_data([sessionpath '\'...
    'all_channels.events']);

% Check the channels that contain TTLs
ChannelID = zeros(1,8);
dataCH    = nan(size(data,1),8);
if any(data==0) == 1
    dataCH(:,1)  = ismember(data,0);
    ChannelID(1) = 1;
end
if any(data==1) == 1
    dataCH(:,2)  = ismember(data,1);
    ChannelID(2) = 2;
end
if any(data==2) == 1
    dataCH(:,3)  = ismember(data,2);
    ChannelID(3) = 3;
end
if any(data==3) == 1
    dataCH(:,4)  = ismember(data,3);
    ChannelID(4) = 4;
end
if any(data==4) == 1
    dataCH(:,5)  = ismember(data,4);
    ChannelID(5) = 5;
end
if any(data==5) == 1
    dataCH(:,6)  = ismember(data,5);
    ChannelID(6) = 6;
end
if any(data==6) == 1
    dataCH(:,7)  = ismember(data,6);
    ChannelID(7) = 7;
end
if any(data==7) == 1
    dataCH(:,8)  = ismember(data,7);
    ChannelID(8) = 8;
end

for ii = find(ChannelID)
    CHdataON  = dataCH(:,ii)==1 & info.eventId==1; % eventId = 1 for onset
    CHdataON  = find(CHdataON);
    CHdataOFF = dataCH(:,ii)==1 & info.eventId==0; % eventId = 0 for offset 
    CHdataOFF = find(CHdataOFF);
    sizedataON = size(CHdataON,1);
    
    % Exclude PulsePal timestamp generated by tagging session (only 1 /
    % tagging session in the beginning of tagging)
    while CHdataON(2) - CHdataON(1) > 10   
        CHdataON(1) = [];               
    end
    if length(CHdataON) < sizedataON - 5    % too many TTLs dropped: suspicious
        error('Error in event convertion.')
    end
    
    CHname = ['CH' num2str(ii)];
    TTL.(CHname) = timestamps(CHdataON);
    pulseon.(CHname)  = timestamps(CHdataON); 
    pulseoff.(CHname) = timestamps(CHdataOFF);  %#ok<FNDSB>
end

sessionpath2 = sessionpath(1:end-26);

% Save bpod TTLs
fnm = fullfile(sessionpath2,'TTLs.mat');
save(fnm,'TTL');

% Save PulsePal TTLs
fnm = fullfile(sessionpath2,'TTLOnOff.mat');
save(fnm,'pulseon','pulseoff');
end