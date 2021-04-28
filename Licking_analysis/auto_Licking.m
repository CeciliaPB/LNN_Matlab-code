
% Read file
path = 'L:\M Kinga\BEHAVIOR\2019 nyár\Cartpt-cre HOMOzigóta 190521\Licking\aborted';
d = dir('*.rhd');
files = {d.name};
n = 600; % time to analize in s, ending at the end of the recording

TabLicks = table;
for ii = 1:length(files)

    file = files{ii};
    read_Intan_RHD2000_file(file,path) % Run and select the .rhd file to open

    lick = board_adc_data(1,:);
    time = t_board_adc;
    maxt = max(time);
    mint = maxt-n;
    time2 = time(time > mint);
    lick2 = lick(time > mint);

    figure('name','Lickport');
    md = median(lick2);
    sd = std(lick2);
    plot(time2, lick2);
    hold on; 
    plot([0 max(time2)],[md md],'LineWidth',1.5);
    plot([0 max(time2)],[md-sd md-sd]);
    plot([0 max(time2)],[md-sd*2 md-sd*2]);
    plot([0 max(time2)],[md-sd*3 md-sd*3]);
    legend('Licking','md','sd','sd*2','sd*3','Location','NorthEastOutside');
    hold off;
    name = file(:,1:end-4);
    saveas(gcf,name,'jpg');
    close gcf;
    
    Th1 = md-sd; 
    Th2 = md-sd*2;
    Th3 = md-sd*3;
    
    numlicks = double(lick2 < Th1);
    licking_time1 = sum(numlicks)/length(time)*100;
    
    numlicks = double(lick2 < Th2);
    licking_time2 = sum(numlicks)/length(time)*100;
    
    numlicks = double(lick2 < Th3);
    licking_time3 = sum(numlicks)/length(time)*100;
    
    TabLicks.file(ii,:) = file;
    TabLicks.sd(ii,:)   = licking_time1;
    TabLicks.sd2(ii,:)  = licking_time2;
    TabLicks.sd3(ii,:)  = licking_time3;

end

writetable(TabLicks,'Results.xls');

%% Lickport 1

% path = cd;
% d = dir('*.rhd');
% files = {d.name};
% n = 600; % time to analize in s, ending at the end of the recording

TabLicks = table;
for ii = 1:size(files,1)

    file = files(ii,:);
    read_Intan_RHD2000_file(file,path) % Run and select the .rhd file to open

    if size(board_adc_data,1) == 1
        lick = board_adc_data(1,:);
    elseif size(board_adc_data,1) == 2
        lick = board_adc_data(2,:);
    end
    time = t_board_adc;
    maxt = max(time);
    mint = maxt-n;
    time2 = time(time > mint);
    lick2 = lick(time > mint);

    clear fit
    fs = frequency_parameters.amplifier_sample_rate;
    bp = [0.001 50]/(fs/2);   %bandpass
    [l,h] = butter(1,bp);
    f = filter(l,h,lick);
    fit = fit(transpose(time),transpose(f),'exp1');
    ff0 = f - fit.a*exp(fit.b*time);
    md = median(ff0);
    sd = std(ff0);
    plot(time, ff0);
    hold on; 
    plot([0 max(time)],[md md],'LineWidth',1.5);
    plot([0 max(time)],[md-sd md-sd]);
    plot([0 max(time)],[md-sd*3 md-sd*3]);
    legend('Licking_F','md','sd','sd*3');
    hold off;
    name = [file(:,1:end-4) '_filt'];
    saveas(gcf,name,'jpg');
    close gcf;
    
    Th1 = md-sd; 
    Th2 = md-sd*2;
    Th3 = md-sd*3;
    
    numlicks = double(ff0 < Th1);
    licking_time1 = sum(numlicks)/length(time)*100;
    
    numlicks = double(ff0 < Th2);
    licking_time2 = sum(numlicks)/length(time)*100;
    
    numlicks = double(ff0 < Th3);
    licking_time3 = sum(numlicks)/length(time)*100;
    
    TabLicks.file(ii,:) = file;
    TabLicks.sd(ii,:)   = licking_time1;
    TabLicks.sd2(ii,:)  = licking_time2;
    TabLicks.sd3(ii,:)  = licking_time3;

end

writetable(TabLicks,'Results_filt.xls');
