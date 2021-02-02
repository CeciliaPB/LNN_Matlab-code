% Pre-assign tables and names ---------------------------------------------
SK = 'Str_SK_MUA';
WN = 'Str_WN_MUA';
PT = 'Str_PT_MUA';

%% Shock data --------------------------------------------------------------
SKtab = load('Neurons_list.mat',SK);
tab     = SKtab.(SK);
bins = 200;

folder = tab.Folder;
Shock = table;
for kk = 1:length(folder)
    
    % Go to folder and load shock TTL
    cd(folder(kk,:));
    
    if kk == 1 
        load('TTLs.mat','TTL');
    elseif min(folder(kk,:) == folder(kk-1,:)) == 1
    elseif min(folder(kk,:) == folder(kk-1,:)) == 0
        clearvars TTL
        load('TTLs.mat','TTL');
    end
    
    neuron = dir(['*' num2str(tab.GR(kk)) '_' num2str(tab.nr(kk)) '.mat']);
    load(neuron.name,'psth0s_shock');
    
    % Shock var to save shock data: FirstSpike and AllSpike
    Shock.Neuron = ['Str_' num2str(kk)];
    Shock.mDelay(kk) = psth0s_shock.FstSpk_means(1,1);
    Shock.mdDelay(kk) = psth0s_shock.FstSpk_means(1,2);
    Shock.stdDelay(kk) = psth0s_shock.FstSpk_means(1,3);
    Shock.numRespNeu(kk) = psth0s_shock.FstSpk_means(1,4)/length(TTL)*100;
    
    [n,edges,bin] = histcounts(psth0s_shock.AllSpk(:,:),bins);
    A = (n == max(n));
    edges(A)
    
    % Plots
    histogram(psth0s_shock.AllSpk(:,:),bins)
    plot(edges(1:bins),n)
end
 


 