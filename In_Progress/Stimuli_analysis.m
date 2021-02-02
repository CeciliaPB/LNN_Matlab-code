
% Pre-assign tables and names ---------------------------------------------
SK = 'BA_SK_MUA';
WN = 'BA_WN_MUA';
PT = 'BA_PT_MUA';

%% Shock data -------------------------------------------------------------
SKtab = load('Neurons_list.mat',SK);
tab   = SKtab.(SK);
bins  = 200;
fig   = 0;

folder  = tab.Folder;
Shock   = table;
for kk = 1:size(folder,1)
    
    % Go to folder and load shock TTL
    cd(folder(kk,:));
    
    if kk == 1 
        load('TTLs.mat','TTL');
    elseif min(folder(kk,:) == folder(kk-1,:)) == 1
    elseif min(folder(kk,:) == folder(kk-1,:)) == 0
        clearvars TTL
        load('TTLs.mat','TTL');
    end
    
    neuron = dir(['GR' num2str(tab.GR(kk)) '_' num2str(tab.nr(kk)) '.mat']);
    if tab.GR(kk) >= 9
        neuron = dir(['PR',num2str(tab.GR(kk)), '_' num2str(tab.nr(kk)) '.mat']);
    end
    variableInfo = who('-file',neuron.name);
    if ismember('psth0s_shock', variableInfo) == 1
        load(neuron.name,'psth0s_shock');
    elseif ismember('psth0_SK', variableInfo) == 1 
        load(neuron.name,'psth0_SK');
        psth0s_shock = psth0_SK;
    elseif ismember('psth0_shock', variableInfo) == 1 
        load(neuron.name,'psth0_shock');
        psth0s_shock = psth0_shock;
    else
        disp(['Loading error in neuron ',SK(1:3) '_' num2str(kk,'%02.0f')])
    end
    
    % Shock var to save shock data: FirstSpike and AllSpike
    Shock.Neuron(kk,:)      = [SK(1:3) '_' num2str(kk,'%02.0f')];
    Shock.mDelay1st(kk)     = psth0s_shock.FstSpk_means(1,1);
    Shock.mdDelay1st(kk)    = psth0s_shock.FstSpk_means(1,2);
    Shock.stdDelay1st(kk)   = psth0s_shock.FstSpk_means(1,3);
    Shock.numRespNeu(kk)    = psth0s_shock.FstSpk_means(1,4)/length(TTL)*100;
    
    [n,edges,~] = histcounts(psth0s_shock.AllSpk(:,:),bins);
        A =(edges>0);
        A = A(1:bins);
        B = n(A);
        C = (B == max(B));
        D = edges(end-length(C)+1:end);
        E = D(C);
        if isempty(E) == 1 
    Shock.mDelayAll(kk)     = NaN;        
        else 
    Shock.mDelayAll(kk)     = E(1);
        end
        A =(edges>0);
        B = edges(A(1:200));
        C = n(A(1:200));
    Shock.stdDelayAll(kk)   = std(B,C);
    
    % Plots
    if fig == 1
        histogram(psth0s_shock.AllSpk(:,:),bins);
        hold on;
        plot(edges(1:bins),n);
        hold off;
    end
end

Shock.MUA  = tab.MUA;
Shock.Resp = tab.Response;

%% White Noise data -------------------------------------------------------
WNtab = load('Neurons_list.mat',WN);
tab   = WNtab.(WN);
bins  = 200;
fig   = 0;

folder  = tab.Folder;
WhiteN  = table;
for kk = 1:size(folder,1)
    
    % Go to folder and load WN TTL
    cd(folder(kk,:));
    
    if kk == 1 
        load('TTLs.mat');
    elseif min(folder(kk,:) == folder(kk-1,:)) == 1
    elseif min(folder(kk,:) == folder(kk-1,:)) == 0
        clearvars -except tab bins fig folder WhiteN Shock kk SK WN PT
        load('TTLs.mat');
    end
    
    % Choose variable to use as stimulus
    if exist('NewWN','var')
        stim = NewWN;
    elseif exist('NewWS','var')
        stim = NewWS;
    elseif exist('TTL_WN','var')
        stim = TTL_WN;
    elseif exist('TTL_WN_video','var')
        stim = TTL_WN_video;
    elseif exist('sound','var')
        stim = sound;
    end
    
    neuron = dir(['GR' num2str(tab.GR(kk)) '_' num2str(tab.nr(kk)) '.mat']);
    if tab.GR(kk) >= 9
        neuron = dir(['PR',num2str(tab.GR(kk)), '_' num2str(tab.nr(kk)) '.mat']);
    end
    variableInfo = who('-file',neuron.name);
    if ismember('psth1s_whiteN', variableInfo) == 1
        load(neuron.name,'psth1s_whiteN');
    elseif ismember('psth1_WN', variableInfo) == 1 
        load(neuron.name,'psth1_WN');
        psth1s_whiteN = psth1_WN;
    else
        disp(['Loading error in neuron',WN(1:3) '_' num2str(kk,'%02.0f')])
    end
    
    % WhiteN var to save shock data: FirstSpike and AllSpike
    WhiteN.Neuron(kk,:)     = [WN(1:3) '_' num2str(kk,'%02.0f')];
    WhiteN.mDelay1st(kk)    = psth1s_whiteN.FstSpk_means(1,1);
    WhiteN.mdDelay1st(kk)   = psth1s_whiteN.FstSpk_means(1,2);
    WhiteN.stdDelay1st(kk)  = psth1s_whiteN.FstSpk_means(1,3);
    WhiteN.numRespNeu(kk)   = psth1s_whiteN.FstSpk_means(1,4)/length(stim)*100;
    
    [n,edges,~] = histcounts(psth1s_whiteN.AllSpk(:,:),bins);
        A =(edges>0);
        A = A(1:bins);
        B = n(A);
        C = (B == max(B));
        D = edges(end-length(C)+1:end);
        E = D(C);
        if isempty(E) == 1 
    WhiteN.mDelayAll(kk)    = NaN;        
        else 
    WhiteN.mDelayAll(kk)    = E(1);
        end
        A =(edges>0);
        B = edges(A(1:200));
        C = n(A(1:200));
    WhiteN.stdDelayAll(kk)  = std(B,C);
    
    % Plots
    if fig == 1
        histogram(psth1s_whiteN.AllSpk(:,:),bins);
        hold on;
        plot(edges(1:bins),n);
        hold off;
    end
end

WhiteN.MUA  = tab.MUA;
WhiteN.Resp = tab.Response;

%% Pure Tone data ---------------------------------------------------------
PTtab = load('Neurons_list.mat',PT);
tab   = PTtab.(PT);
bins  = 200;
fig   = 0;

folder  = tab.Folder;
PureT   = table;
for kk = 1:size(folder,1)
    
    % Go to folder and load shock TTL
    cd(folder(kk,:));
    
    if kk == 1 
        load('TTLs.mat');
    elseif min(folder(kk,:) == folder(kk-1,:)) == 1
    elseif min(folder(kk,:) == folder(kk-1,:)) == 0
        clearvars -except tab bins fig folder PureT WhiteN Shock kk SK WN PT
        load('TTLs.mat');
    end
    
    % Choose variable to use as stimulus
    if exist('NewPT','var')
        stim = NewPT;
    elseif exist('TTL_PT','var')
        stim = TTL_PT;
    elseif exist('TTL_PT_video','var')
        stim = TTL_PT_video;
    elseif exist('sound','var')
        stim = sound;
    end
    
    neuron = dir(['GR' num2str(tab.GR(kk)) '_' num2str(tab.nr(kk)) '.mat']);
    if tab.GR(kk) >= 9
        neuron = dir(['PR',num2str(tab.GR(kk)), '_' num2str(tab.nr(kk)) '.mat']);
    end
    variableInfo = who('-file',neuron.name);
    if ismember('psth1s_pureT', variableInfo) == 1
        load(neuron.name,'psth1s_pureT');
    elseif ismember('psth1_PT', variableInfo) == 1 
        load(neuron.name,'psth1_PT');
        psth1s_pureT = psth1_PT;
    else
        disp(['Loading error in neuron',PT(1:3) '_' num2str(kk,'%02.0f')])
    end
    
    % Shock var to save shock data: FirstSpike and AllSpike
    PureT.Neuron(kk,:)      = [PT(1:3) '_' num2str(kk,'%02.0f')];
    PureT.mDelay1st(kk)     = psth1s_pureT.FstSpk_means(1,1);
    PureT.mdDelay1st(kk)    = psth1s_pureT.FstSpk_means(1,2);
    PureT.stdDelay1st(kk)   = psth1s_pureT.FstSpk_means(1,3);
    PureT.numRespNeu(kk)    = psth1s_pureT.FstSpk_means(1,4)/length(stim)*100;
    
    [n,edges,~] = histcounts(psth1s_pureT.AllSpk(:,:),bins);
        A =(edges>0);
        A = A(1:bins);
        B = n(A);
        C = (B == max(B));
        D = edges(end-length(C)+1:end);
        E = D(C);
        if isempty(E) == 1 
    PureT.mDelayAll(kk)     = NaN;        
        else 
    PureT.mDelayAll(kk)     = E(1);
        end
        A =(edges>0);
        B = edges(A(1:200));
        C = n(A(1:200));
    PureT.stdDelayAll(kk)   = std(B,C);
    
    % Plots
    if fig == 1
        histogram(psth1s_pureT.AllSpk(:,:),bins);
        hold on;
        plot(edges(1:bins),n);
        hold off;
    end
end

PureT.MUA  = tab.MUA;
PureT.Resp = tab.Response;

%% Plots ------------------------------------------------------------------

% Location = strfind(Neurons_List.Location,'LA');
% SK = Neurons_List.Shock == -1 ;
% for ii = 1:length(Location)
%     try
%         A(ii) = cell2mat(Location(ii));
%     catch
%         A(ii) = 0;
%     end
% end
% LAsk = sum(SKdata == 1 & transpose(Location == 1));

% SKdata = Shock.mdDelay1st(Shock.Resp == 1 & Shock.MUA == 0); 
% WNdata = WhiteN.mdDelay1st(WhiteN.Resp == 1 & WhiteN.MUA == 0);
% PTdata = PureT.mdDelay1st(PureT.Resp == 1 & PureT.MUA == 0);
SKdata = Shock.mdDelay1st(Shock.Resp == 1 & Shock.MUA == 1); 
WNdata = WhiteN.mdDelay1st(WhiteN.Resp == 1 & WhiteN.MUA == 1);
PTdata = PureT.mdDelay1st(PureT.Resp == 1 & PureT.MUA == 1);

% Shock: Color 255 204 51
y = 0.9 + (1.1-0.9).*rand(length(SKdata),1); 
SKcolor = 1/255*[255,204,51];
SKscat = 1/255*[102,102,102];
hold on; boxplot(SKdata*1000,'symbol','','Color',SKcolor,'Labels',{'Shock'},'Widths',0.2)
h = findobj(gca,'Tag','Box');
uistack(patch(get(h,'XData'),get(h,'YData'),SKcolor,'FaceAlpha',.5,'EdgeColor',SKcolor,'LineWidth',2),'bottom');
set(findobj(gcf,'-regexp','Tag','\w*Whisker'),'LineStyle','-','LineWidth',2)
set(findobj(gcf,'-regexp','Tag','\w*Adjacent'),'LineWidth',2)
set(findobj(gcf,'-regexp','Tag','Median'),'LineWidth',3,'Color','r');
scatter(y,SKdata*1000,'MarkerFaceColor',SKscat,'MarkerEdgeColor',SKscat,'LineWidth',1)
ylim([-5 300]);
ylabel('Time (ms)');
set(gca,'fontname','arial');
mdSKdata = nanmedian(SKdata);
sdSKdata = nanstd(SKdata);
text(0.55,250,['Delay: ' num2str(mdSKdata*1000) '\pm' num2str(sdSKdata*1000)]);
text(0.55,270,['N = ' num2str(size(SKdata,1))]);

% White Noise: Color 0 153 204
figure;
y = 0.9 + (1.1-0.9).*rand(length(WNdata),1);
WNcolor = 1/255*[0, 51, 204];
SKscat = 1/255*[102,102,102];
hold on; boxplot(WNdata*1000,'symbol','','Color',WNcolor,'Labels',{'White noise'},'Widths',0.2)
h = findobj(gca,'Tag','Box');
uistack(patch(get(h,'XData'),get(h,'YData'),WNcolor,'FaceAlpha',.5,'EdgeColor',WNcolor,'LineWidth',2),'bottom');
set(findobj(gcf,'-regexp','Tag','\w*Whisker'),'LineStyle','-','LineWidth',2)
set(findobj(gcf,'-regexp','Tag','\w*Adjacent'),'LineWidth',2)
set(findobj(gcf,'-regexp','Tag','Median'),'LineWidth',3,'Color','r');
scatter(y,WNdata*1000,'MarkerFaceColor',SKscat,'MarkerEdgeColor',SKscat,'LineWidth',1)
ylim([-5 300]);
ylabel('Time (ms)');
set(gca,'fontname','arial');
mdWNdata = nanmedian(WNdata);
sdWNdata = nanstd(WNdata);
text(0.55,250,['Delay: ' num2str(mdWNdata*1000) '\pm' num2str(sdWNdata*1000)]);
text(0.55,270,['N = ' num2str(size(WNdata,1))]);

% Pure Tone: Color 0 153 204
figure;
y = 0.9 + (1.1-0.9).*rand(length(PTdata),1);
PTcolor = 1/255*[0, 153, 204];
SKscat = 1/255*[102,102,102];
hold on; boxplot(PTdata*1000,'symbol','','Color',PTcolor,'Labels',{'Pure tone'},'Widths',0.2)
h = findobj(gca,'Tag','Box');
uistack(patch(get(h,'XData'),get(h,'YData'),PTcolor,'FaceAlpha',.5,'EdgeColor',PTcolor,'LineWidth',2),'bottom');
set(findobj(gcf,'-regexp','Tag','\w*Whisker'),'LineStyle','-','LineWidth',2)
set(findobj(gcf,'-regexp','Tag','\w*Adjacent'),'LineWidth',2)
set(findobj(gcf,'-regexp','Tag','Median'),'LineWidth',3,'Color','r');
scatter(y,PTdata*1000,'MarkerFaceColor',SKscat,'MarkerEdgeColor',SKscat,'LineWidth',1)
ylim([-5 300]);
ylabel('Time (ms)');
set(gca,'fontname','arial');
mdPTdata = nanmedian(PTdata);
sdPTdata = nanstd(PTdata);
text(0.55,250,['Delay: ' num2str(mdPTdata*1000) '\pm' num2str(sdPTdata*1000)]);
text(0.55,270,['N = ' num2str(size(PTdata,1))]);
