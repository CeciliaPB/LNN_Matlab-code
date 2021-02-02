NeuronsResp = table;
SK = abs(NeuronsResp.SK).*transpose(t);
WN = abs(NeuronsResp.WN).*transpose(t);
PT = abs(NeuronsResp.PT).*transpose(t);

%% Analysing shock

Nfolder = Str_SK_MUA.Folder;

% All shock responsive neurons --------------------------------------------
Shock = table;
for ii = 1:length(Nfolder)
    cd(Nfolder(ii,:));
    ID = [num2str(Str_SK_MUA.GR(ii)), '_', num2str(Str_SK_MUA.nr(ii))];
    
    try 
        load(['GR',ID,'.mat'], 'psth0s_shock');
    catch
        load(['TT',ID,'.mat'], 'psth0s_shock');
    end
 
    Shock.Neuron = sprintf('Str_%02d',ii);
    Shock.mDelay1st(ii) = psth0s.FstSpk_means(1,1);
    Shock.mdDelay1st(ii) = psth0s.FstSpk_means(1,2);
    Shock.stdDelay1st(ii) = psth0s.FstSpk_means(1,3);
    Shock.numRespNeu(ii) = psth0s.FstSpk_means(1,4)/length(TTL)*100;
end


y = ones(length(SK(:,1)),1);
hold on; boxplot(SK(:,1),'Labels',{'Shock'},'Widths',0.05)
scatter(y,SK(:,1),'filled','MarkerEdgeColor','flat','LineWidth',2)
ylim([0.0 0.18]);
figure;
y = ones(length(WN(:,1)),1);
hold on; boxplot(WN(:,1),'Labels',{'White noise'},'Widths',0.05)
scatter(y,WN(:,1),'filled','MarkerEdgeColor','flat','LineWidth',2)
ylim([0.0 0.18]);
figure;
y = ones(length(PT(:,1)),1);
hold on; boxplot(PT(:,1),'Labels',{'Pure tone'},'Widths',0.05)
scatter(y,PT(:,1),'filled','MarkerEdgeColor','flat','LineWidth',2)
ylim([0.0 0.18]);

WN1 = WN(:,1)<0.1;
WN2 = WN(:,1)>0.1;
WN1 = WN(WN1);
WN2 = WN(WN2);
y = ones(length(WN1),1);
hold on; boxplot(WN1,'Labels',{'Pure tone'},'Widths',0.05)
scatter(y,WN1,'filled','MarkerEdgeColor','flat','LineWidth',2)
ylim([0.0 0.18]);
figure;
y = ones(length(WN2),1);
hold on; boxplot(WN2,'Labels',{'Pure tone'},'Widths',0.05)
scatter(y,WN2,'filled','MarkerEdgeColor','flat','LineWidth',2)
ylim([0.0 0.18]);
