function cb_mksessionfolder(AnimalID,Task)
%   cb_mksessionfolder makes the folders in CellBase format for easy use in
% 	CB.
%   INPUTS:
%       - AnimalID: string, code of the animal. Ex. 'BAd01'.
%       - Task: string, behavioural protocol. Ex. 'Gambling'.
%   It relies on CellBase data handling system.
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2021
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

B = dir([AnimalID '_' file(end-18:end-15) '-' file(end-14:end-13)...
    '-' file(end-12:end-11) '*']);
names2 = {B.name}';
B = cell2mat({B.isdir}');
for kk = 1:length(names2)
    for ii = 1:length(names)
        file = names{ii,1};
        folder = [file(20:end-11) 'a'];
        folder2 = [file(20:end-11) 'b'];
        if kk == 1 && B(kk) == 1 && exist('folder','var')
            movefile(names2{kk,1}, folder);
        elseif kk > 1 && B(kk) == 1 && exist('folder2','var')
            movefile(names2{kk,1}, folder2);
        else
            disp('No recording found.');
        end
    end
end

end
