%% Calculating firing rate
% [~, time, ~] = load_open_ephys_data('100_CH1.continuous');
% TT = ((TS(:,1)/10000)-(min(time))+Dt);

iFRn = [];
for ii = 1:27
    time = [NewPT(ii)-0.5,NewPT(ii)+0.5];
    timestamps = TT(TT > time(1) & TT < time(2));
    [iFR,~] = instFiringRate(timestamps,time,100);
    iFRn(ii,:) = iFR;
end
t = linspace(-0.5,0.5,size(iFRn,2));
figure;plot(t,iFRn)
iFRm = mean(iFRn);
hold on;
plot(t,iFRm,'k','Linewidth',2)
plot([0 0],[0 600],'r','Linewidth',2)

%% Extract data from wf

Tab    = Lick;
try
    folder = Tab.folder;
catch
    folder = Tab.Folder;
end
Data   = table;

for kk = 1:size(Tab.folder,1)
    
    % Go to folder and load shock TTL
    cd(folder(kk,:));
    
    % Find neuron
    NeuronID = dir(['*',num2str(Tab.GR(kk)),'_',num2str(Tab.nr(kk)),'.mat']);
    NeuronID = NeuronID.name;
    
%     load(NeuronID,'wf');
%     
%     if exist('wf','var') == 0
        % calculate waveforms
        A = dir(['*',num2str(Tab.GR(kk)),'.mat']);
        B = max(cell2mat({A.bytes})); % Select the largest (will have TimeStamps and WaveForms) 
        C = find(cell2mat({A.bytes}) == B);
        D = {A.name};
        E = D(C);
        GroupID = E{1};
        load(GroupID,'TimeStamps','WaveForms');
        
        wf = waveform_analysis(GroupID,NeuronID,'features');
        
        save(NeuronID,'wf','-append');
%     end
    
    % load the variables inside wf
    names = fieldnames(wf);
    for ii=1:length(names)
        eval([names{ii} '=wf.' names{ii} ]);
    end
    
    if pk1(1,2)>pk2(1,2) % In case there are no two valleys
        b   = a;
        a   = 0;
        pk3 = pk1;
        pk1 = [0,0];
        d   = abs(c);
        c   = NaN;
        f   = e;
        e   = NaN;
        first = NaN;
    end    
    
    % Calculations
    sym     = a/(a+b);
    amp1    = abs(e); % pk2-pk1
    amp2    = abs(f); % pk3-pk2
    pkdist1 = (15 * 0.5 )/c; % in ms, distance pk1 to pk2
    pkdist2 = (15 * 0.5 )/d; % in ms, distance pk2 to pk3
    halfwidth = (15 * 0.5 )/g; % in ms, halfwith pk2
    pkwidth1  = (15 * 0.5 )/first; % in ms, width pk1
    pkwidth2  = (15 * 0.5 )/last; % in ms, width pk3
    
    Data.neuron(kk,:) = ([folder(kk,13:end), '_',NeuronID(1,end-8:end-4)]);
    Data.sym(kk)      = sym;
    Data.amp1(kk)     = amp1;
    Data.amp2(kk)     = amp2;
    Data.pkdist1(kk)  = pkdist1;
    Data.pkdist2(kk)  = pkdist2;
    Data.halfwidth(kk) = halfwidth;
    Data.pkwidth1(kk)  = pkwidth1;
    Data.pkwidth2(kk)  = pkwidth2;
    
    clearvars -except Data Tab folder Exploration* Groom* Lick* Rest* Str_PT_MUA Str_WN_MUA PT* WN*
end

%% Statistical analysis
% amp2

Var = 'dist';
KWData = nan(29,6);
KWData(1:length(ExplorationData.(Var)),1) = [ExplorationData.(Var)];
KWData(1:length(RestData.(Var)),2)        = [RestData.(Var)];
KWData(1:length(GroomData.(Var)),3)       = [GroomData.(Var)];
KWData(1:length(LickData.(Var)),4)        = [LickData.(Var)];
KWData(1:length(WNData.(Var)),5)          = [WNData.(Var)];
KWData(1:length(PTData.(Var)),6)          = [PTData.(Var)];

% % Normality h = 0, normal; h = 1, not normal.
% for ii = 1:6
% h(ii) = kstest(KWData(:,ii));
% end

% % Homoscedasticity p < 0.05, not same variance; p > 0.05, same variance.
% pvar = vartestn(KWData(:,[1:3,5:6]));

[p,tbl,stats] = kruskalwallis(KWData(:,5:6));
% c = multcompare(stats);
% [p,tbl,stats] = kruskalwallis((KWData(:,[1:3,5:6])).^2);
% c = multcompare(stats);
% A = zero2nan(nansum(KWData(:,1:3),2));
% B = zero2nan(nansum(KWData(:,5:6),2));
% p2 = kruskalwallis(([A,B]).^2);

%% Plots
PlotData = KWData;

% Pure Tone: Color 0 153 204
figure;
y(:,1) = 0.9 + (1.1-0.9).*rand(length(PlotData),1);
y(:,2) = 1.9 + (2.1-1.9).*rand(length(PlotData),1);
y(:,3) = 2.9 + (3.1-2.9).*rand(length(PlotData),1);
y(:,4) = 3.9 + (4.1-3.9).*rand(length(PlotData),1);
y(:,5) = 4.9 + (5.1-4.9).*rand(length(PlotData),1);
y(:,6) = 5.9 + (6.1-5.9).*rand(length(PlotData),1);
Color = 1/255*[0, 153, 204];
SKscat = 1/255*[102,102,102];
hold on; boxplot(PlotData,'symbol','','Color',Color,'Labels',{'Exploring','Resting','Grooming','Licking','White Noise','Pure Tone'},'Widths',0.5)
h = findobj(gca,'Tag','Box');
% uistack(patch(get(h,'XData'),get(h,'YData'),Color,'FaceAlpha',.5,'EdgeColor',Color,'LineWidth',2),'bottom');
set(findobj(gcf,'-regexp','Tag','\w*Whisker'),'LineStyle','-','LineWidth',2)
set(findobj(gcf,'-regexp','Tag','\w*Adjacent'),'LineWidth',2)
set(findobj(gcf,'-regexp','Tag','Median'),'LineWidth',3,'Color','r');
plot(y,PlotData,'o','MarkerFaceColor',SKscat,'MarkerEdgeColor',SKscat,'LineWidth',1)
ylim([-0.1 1]);
ylabel('Time (ms)');
set(gca,'fontname','arial');

