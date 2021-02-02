
% Read file
RHD_list = dir(['*PrL2' '*.rhd']);
RHD_list = {RHD_list.name}';

%%
for kk = 1:10
name  = RHD_list{kk};
read_Intan_RHD2000_file(name,cd)

% read_Intan_RHD2000_file % Run and select the .rhd file to open

lick_0 = board_adc_data(1,:);
time = t_board_adc;

% Lickport 0

figure('name','Lickport 0');
% Subplot 1: Raw data with median and sd values
md = median(lick_0);
sd = std(lick_0);
% subplot(2,1,1); 
plot(time, lick_0);
hold on; 
plot([0 max(time)],[md md],'LineWidth',1.5);
plot([0 max(time)],[md+sd md+sd]);
plot([0 max(time)],[md-sd md-sd]);
plot([0 max(time)],[md+sd*2 md+sd*2]);
legend('Licking','md','sd','-sd','sd*2');
hold off;
% % Subplot 2: Filtered data, for a better detection
% clear fit
% bp = [0.001 50]/(frequency_parameters.amplifier_sample_rate/2);   %bandpass
% [l,h] = butter(1,bp);
% f = filter(l,h,lick_0);
% fit=fit(transpose(time),transpose(f),'exp1');
% ff0 = f - fit.a*exp(fit.b*time);
% mdf = median(ff0);
% sdf = std(ff0);
% subplot(2,1,2); plot(time, ff0);
% hold on; 
% plot([0 max(time)],[mdf mdf],'LineWidth',1.5);
% plot([0 max(time)],[mdf-sdf mdf-sdf]);
% plot([0 max(time)],[mdf+sdf mdf+sdf]);
% legend('Licking_F','mdf','sdf');

% Thresholding
Th = input('Threshold value? ENTER VALUE:'); % You can enter numeric value or refer to an existing variable
% Sg = input('Signal? (raw or filt) ENTER VAR NAME: ','s'); % Enter name of variable
% if Sg == 'raw'
%     Signal = lick_0;
% elseif Sg ~= 'raw'
%     Signal = ff0;
% 
% end
Signal = lick_0;

if Th < md
licking_inout = Signal < Th;
elseif Th > md
licking_inout = Signal > Th;
end
licking_inout = double(licking_inout);
% plot(time,licktime_0*0.1); 
hold on;
plot(time,licking_inout*0.4);
hold off;

% lick_time_100 = sum(licking_inout)/length(time)*100; 

idx = find(licking_inout == 1);
inout2 = licking_inout;
inout2(end) = 0;
n = ceil(0.0002*20000); % 0.0002 time duration of events (licks)

for ii = 1:length(idx)-n*3
    if inout2(idx(ii)) == 0

    elseif inout2(idx(ii)) == 1
        A = sum(inout2(idx(ii):idx(ii)+n*3));
        if A >= ceil(n*2.5)
            inout2(idx(ii)) = 1;
        elseif A ~= ceil(n*2.5)
            inout2(idx(ii)) = 0;
        end
    end
end

hold on;
plot(time,inout2*0.3);
legend('Licking','md','sd','-sd','sd*2','Licks detected','Actual Licks')
title([name(1:8),' Licking'])
hold off;

licking_inout = inout2;
lick_time_100 = sum(licking_inout)/length(time)*100;

licking_inout2 = zeros(1,length(licking_inout));
   for ii = 1:length(licking_inout)-2000
    A = sum(licking_inout(ii:ii+2000));
    if A > 500
      licking_inout2(ii) = 1;
    elseif A <= 500
      licking_inout2(ii) = 0;
    end
   end
licking_inout2 = downsample(licking_inout2,1000);
disp('DONE!');

if kk == 12 || kk == 15
    save([name(1:7),'_Test_Track.mat'],'lick_time_100','licking_inout','licking_inout2','-append');
else
    save([name(1:8),'_track.mat'],'lick_time_100','licking_inout','licking_inout2','-append');
end

end