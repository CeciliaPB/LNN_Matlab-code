% Exploration neurons
Beh = Sound;

Data = table;
for ii = 1:size(Beh.folder,1)
    
folder = Beh.folder(ii,:);
GR = Beh.GR(ii);
nr = Beh.nr(ii);
neuron = ([folder, filesep, 'GR',num2str(GR),'_',num2str(nr),'.mat']);

TT = load(neuron,'wf');
WF = TT.wf;

names = fieldnames(WF);
for kk=1:length(names)
eval([names{kk} '=WF.' names{kk} ]);
end

if pk1(1,2)>pk2(1,2)
    b = a;
    a = NaN;
    pk3 = pk1;
    pk1 = [0,0];
    d = abs(c);
    c = NaN;
    f = e;
    e = NaN;
    first = NaN;
end    

sym = a/(a+b);
amp1 = abs(e); % pk2-pk1
amp2 = abs(f); % pk3-pk2
pkdist1 = (15 * 0.5 )/c; % in ms, distance pk1 to pk2
pkdist2 = (15 * 0.5 )/d; % in ms, distance pk2 to pk3
halfwidth = (15 * 0.5 )/g; % in ms, halfwith pk2
pkwidth1 = (15 * 0.5 )/first; % in ms, width pk1
pkwidth2 = (15 * 0.5 )/last; % in ms, width pk3

A = folder(1,end-4:end);

Data.neuron(ii,:) = ([A, '_TT',num2str(GR),'_',num2str(nr)]);
Data.sym(ii)  = sym;
Data.amp1(ii) = amp1;
Data.amp2(ii) = amp2;
Data.pkdist1(ii) = pkdist1;
Data.pkdist2(ii) = pkdist2;
Data.halfwidth(ii) = halfwidth;
Data.pkwidth1(ii)  = pkwidth1;
Data.pkwidth2(ii)  = pkwidth2;

end
