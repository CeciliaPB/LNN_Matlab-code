function cb_mksessionfolder(AnimalID,Task)
%   cb_mksessionfolder makes the folders in CellBase format for easy use in
% 	CB.
%   INPUTS:
%       - AnimalID: string, code of the animal. Ex. 'BAd01'.
%       - Task: string, behavioural protocol. Ex. 'Gambling'.
%   It relies on CellBase data handling system.
% -------------------------------------------------------------------------
% Cec�lia Pardo-Bellver, 2021
% Laboratory of Network Neurophysiology
% Institute of Experimantal Medicine, Hungary.
% -------------------------------------------------------------------------

A = dir([AnimalID '_' Task '*.mat']);
names = {A.name}';

for ii = 1:length(names)
    file = names{ii,1};
    folder = [file(20:end-11) 'a'];
    
    if ~exist(folder,'dir') == 1
        mkdir(folder);
        movefile(file, folder);
    elseif ~exist(folder,'dir') == 0
        folder2 = [folder(1:end-1) 'b'];
        mkdir(folder2);
        movefile(file, folder2);
    end
end
end
