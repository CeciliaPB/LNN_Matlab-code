%% From Lickport data

fs = 30000;
lick = data(2,:);

figure('name','Lickport 0');
% Subplot 1: Raw data with median and sd values
clear fit
bp = [0.001 50]/(fs/2);   %bandpass
[l,h] = butter(1,bp);
f = filter(l,h,lick);
fit = fit(transpose(time),transpose(f),'exp1');
ff0 = f - fit.a*exp(fit.b*time);
mdf = median(ff0);
sdf = std(ff0);
plot(time, ff0);
hold on; 
plot([0 max(time)],[mdf mdf],'LineWidth',1.5);
plot([0 max(time)],[mdf-sdf mdf-sdf]);
plot([0 max(time)],[mdf-sdf*3 mdf-sdf*3]);
legend('Licking_F','mdf','sdf','sdf*3');

% Thresholding to detect below Th values (licking)
Th = input('Threshold value? ENTER VALUE:'); 
% Th: You can enter numeric value or refer to an existing variable.
Signal = ff0;

licktime = double(Signal < Th); % Vector, 0 =  above Th; 1 = below Th.
licktime(1,1) = 0;
% plot(time,licktime_0*0.1); % Non-filtered 
plot(time,(licktime-2)*0.2); % Filtered data

% Percentage of time spent licking or trying
lick_time = sum(licktime); 
lick_100 = sum(licktime)/length(time)*100;

% Number of lickings
A = ischange(licktime,'MaxNumChanges',length(time));
B = time(A); 
if mod(length(B),2)==0 
elseif mod(length(B),2)~=0
    B = [B,t2(end)];
end
lick_start = transpose(B(:,1:2:end));  % odd matrix, start lick
lick_end = transpose(B(:,2:2:end)); % even matrix, end lick

NumLicks = length(lick_start);
fprintf(['Number of licks = ', num2str(NumLicks)]);

C = lick_end - lick_start;
D = (C > 0.01);
lick_start(D == 0) = NaN;
lick_end(D == 0) = NaN;
E = abs(double(isnan(lick_start)) - 1);
licks(:,1) = lick_start(logical(E));
E = abs(double(isnan(lick_end)) - 1);
licks(:,2) = lick_end(logical(E));

%% Lick frequency

Lick_freq = [];
ii = 2;
t = lick_start(ii-1,1);

while ii < length(lick_start)
    
    int = lick_start(ii,1) - lick_end(ii-1,1);
    
    if int < 0.5
        t = [t; lick_start(ii,1)];
    elseif int >= 0.5
        freq = length(t)/(t(end) - t (1));
        Lick_freq = [Lick_freq, freq];
        clearvars t
        t = lick_start(ii,1);
    end
    
    ii = ii + 1;
end

if ii == length(lick_start)
    freq = length(t)/(t(end) - t (1));
    Lick_freq = [Lick_freq, freq];
end

%% From the video
tmp = dir('*rackL.mat');
dat = {tmp.name}';
name(1:9,1) = dat(2:10,1);
name(10,1) = dat(1,1);

tmp = dir('*rack.mat');
dat = {tmp.name}';
name2(1:9,1) = dat(2:10,1);
name2(10,1) = dat(1,1);

allLicks = [];
for kk = 1:length(name)

load(name{kk})
load(name2{kk})
lick = trackL(:,1);
time = tL;

L = abs(double(isnan(lick)) - 1);

% As small body entries around the port can be detected, this sets a min
% num of frames (3) to be considered a lick.
n = 3; % num of frames
lick_in = L;
for ii = 1:length(L)-n
    A = nansum(L(ii+1:ii+n,1));
    B = nansum(L(ii+1:ii-n,1));
    if L(ii) == 1 && A >= 1
        lick_in(ii,1) = 1;
    elseif L(ii) == 1 && A == 0
        lick_in(ii,1) = 0;
    elseif L(ii) == 0 && A >= 1 && B >= 1
        lick_in(ii,1) = 1;
    elseif L(ii) == 0 && A >= 2 && B == 0
        lick_in(ii,1) = 0;
    elseif L(ii) == 0 && A == 0
        lick_in(ii,1) = 0;
    end
end

% Number of approaches 
Licks = zeros(length(trigger),1);
for jj = 1:length(trigger)
    fin = trigger(jj);
    
    trDur = (tL>=(fin-150) & tL<=fin); trDur(1,end) = 0;
    trLick = transpose(lick_in) & trDur;
    
    A = ischange(double(trLick),'MaxNumChanges',length(tL));
    B = time(A); 
    if mod(length(B),2)==0 
    elseif mod(length(B),2)~=0
        B = [B,tL(end)];
    end
    lick_start = transpose(B(:,1:2:end));  % odd matrix, start lick
    lick_end = transpose(B(:,2:2:end)); % even matrix, end lick

    NumLicks = length(lick_start);
    
    Licks(jj,1) = NumLicks;
    Licks(jj,2) = sum(lick_end - lick_start);
end

allLicks = [allLicks; Licks];
save(name{kk},'Licks','-append');

end

%% Add to existing table
tmp = dir('*rackL.mat');
dat = {tmp.name}';
C(1:9,1) = dat(2:10,1);
C(10,1) = dat(1,1);

for ii = 1:9:length(dat)*9
   n = ceil(ii/9);
   load(C{n});
   A8.licks(ii:ii+8,1) = Licks(:,1);
   A8.licks_dur(ii:ii+8,1) = Licks(:,2);  
end

% licks = zeros(15,1);
% for ii = 1:3:90
%     B = nanmean(A2.licks(ii:ii+2,:));
%     n = ceil(ii/3);
%     licks(n,1) = B;
% end
% 
% plot(licks(:,1),'-o','LineWidth',2,'MarkerFaceColor','auto');
