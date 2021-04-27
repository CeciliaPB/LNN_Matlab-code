function channel_map(CHfolder, processor, probe)

% This function reads .continuous files from OpenEphys, renames and saves
% them in a new folder (CHfolder).
% 
% INPUTS: 
%   - CHfolder: Name of the folder where the reordered channels will be
%   saved. 
%   - processor: processor number. Default 101.
%   - probe: probe type, the script reorders the channel numbers according
%   to it.
%
% Examples: 
%   channel_map('ReorderCH','110','A32');
% 
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2019
% Laboratory of Network Neurophysiology
% Institute of Experimantal Medicine, Hungary.
%
% MATLAB toolbox: Signal Processing Toolbox.
% -------------------------------------------------------------------------

% Channel order -----------------------------------------------------------
switch probe
    case 'A32'
        channel = [1 17 16 32 3 19 14 30 9 25 10 20 8 24 2 29 7 26 15 21 11 23 12 28 6 18 13 22 5 27 4 31];
    case 'A4x8' % Female view probe/ female view adaptor
        channel = [28 18 23 22 21 27 26 31 30 25 19 20 32 24 17 29 9 14 10 3 8 16 2 1 6 12 13 11 5 15 4 7];
    otherwise
        disp('Unknown probe');
end

% Making and setting folders ----------------------------------------------
if ~isfolder(CHfolder)
    mkdir(CHfolder)
end
currentfolder = pwd;
newfolder = [currentfolder filesep CHfolder];

% This reorders the files so after 1 comes 2 and not 10. Thank you MATLAB.
CHfiles = dir([num2str(processor) '_CH' '*.continuous']);
CHfiles = {CHfiles.name};
CHfiles = natsortfiles(CHfiles);

for ii = 1:length(channel)
   
    ch = channel(ii);
    newfile = [newfolder filesep '200_CH' num2str(ii) '.continuous'];
    oldfile = [currentfolder filesep CHfiles{ch}];
    copyfile(oldfile,newfile);  
end

end

%%
% AnimalID = 'HB6';
% A = dir([AnimalID '_G' '*.mat']);
% names = {A.name}';
% for ii = 1:length(names)
% file = names{ii,1};
% folder = [file(18:end-11) 'a'];
% if ~exist(folder,'dir') == 1
% mkdir(folder);
% movefile(file, folder);
% elseif ~exist(folder,'dir') == 0
% folder2 = [folder(1:end-1) 'b'];
% mkdir(folder2);
% movefile(file, folder2);
% end
% end
 