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

idlist      = {'a' 'b' 'c' 'd' 'e'};
for ii = 1:length(names)
    sessionlist = dir('20*');
    sessionlist = {sessionlist.name}';
    file = names{ii,1};
    session = file(end-18:end-11);
    isSession = contains(sessionlist, session);
    if any(isSession) == 0
        folder = [file(end-18:end-11) idlist{1}];
        mkdir(folder);
        movefile(file, folder);
    elseif any(isSession) == 1
        lastSession = sessionlist(isSession == 1);
        lastSession = lastSession{end};
        lastId = find(contains(idlist, lastSession(end)));
        folder = [file(end-18:end-11) idlist{lastId+1}];
        mkdir(folder);
        movefile(file, folder);       
    end
end

B = dir([AnimalID '_20*']);
names2 = {B.name}';

sessionlist = dir('20*');
sessionlist = {sessionlist.name}';
for kk = 1:length(names2)
    Part1   = names2{kk}(end-18:end-15);
    Part2   = names2{kk}(end-13:end-12);
    Part3   = names2{kk}(end-10:end-9);
    recSession = [Part1, Part2, Part3];
    isSession  = contains(sessionlist, recSession);
    lastRec    = sessionlist(isSession == 1);
    isRec = dir([lastRec{1} '\' AnimalID '_20*']);
    if isempty(isRec) == 0
        disp([lastRec{1} ' already contains a recording.']);
    elseif isempty(isRec) == 1
        folder = lastRec{1};
        movefile(names2{kk,1}, folder);
    else
        disp('No recording found.');
    end
end

end
