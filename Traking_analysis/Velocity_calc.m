%% 1_ load the .mat file and select the video

video_list = dir(['E','*.mp4']);
video_list = {video_list.name}';
ii = 4;
name = video_list{ii};
load([name(1:10),'.mat']);

%% 2_ in case you forgot to save the inArea, calculate it

if ~exist('inArea', 'var') 
    video = name;
    videoSource = vision.VideoFileReader(video); % Read the video
    frameSize = size(step(videoSource));
    videoPlayer = vision.VideoPlayer('Position',[100 100 [frameSize(2),...
    frameSize(1)]+30]); 

    [X,Y] = meshgrid(1:frameSize(2),1:frameSize(1));
    frame = step(videoSource);
    figure('name','Area'); imshow(step(videoSource));
    [x,y] = ginput;
    inArea = inpolygon(X(:),Y(:),x,y);
    inArea = double(reshape(inArea,[frameSize(1) frameSize(2)]));
    frame(~inArea) = 0;
    hold on; imshow(frame); hold off;
end

%% 3_ plot the frame and find the coordinates for the corners

imshow(inArea)
hold on;
plot(track(:,1),track(:,2),'r','LineWidth',1);
% Select the DataTips tool and look for the XY info of the corners, to save
% the plot: saveas(gcf, name(1:end-4), 'jpg');

%% 4_ write the coordinates for the corners

% If you calculate the inArea, the coordinates come from x,y values
% Corners: B ---- C
%          |      |
%          |      |
%          A ---- D
%  A = [33 402];
%  B = [23 3];
%  C = [560 2];
%  D = [570 394]; 

%% 5_ calculate distances (vectors)
% X1 = BC;
% X2 = AD;
% Y1 = AB;
% Y2 = CD;
X1 = sqrt(((B(1)-C(1))^2)+((B(2)-C(2))^2));
X2 = sqrt(((A(1)-D(1))^2)+((A(2)-D(2))^2));
Y1 = sqrt(((A(1)-B(1))^2)+((A(2)-B(2))^2));
Y2 = sqrt(((C(1)-D(1))^2)+((C(2)-D(2))^2));

X = mean([X1,X2]);
Y = mean([Y1,Y2]);

%% 6_ calculate new track values (change px to cm)

trackX = (track(:,1) * 24.5)/X; % Where 24.5 is the measured X distance in cm
trackY = (track(:,2) * 19.4)/Y; % Where 19.4 is the measured Y distance in cm

%% 7_ calculate the time vector

info = audioinfo(name);
time = linspace(0,info.Duration,length(track));
time = transpose(time); 

%% 8_ Calculate Speed, set the point distances between frames

Distance = zeros(length(trackX)-1,1);
for ii = 1:length(Distance)
frame1 = [trackX(ii,1) trackY(ii,1)]; % Extracts xy for the frame
frame2 = [trackX(ii+1,1) trackY(ii+1,1)]; % Extracts xy for the NEXT frame

Dist_vect = sqrt(((frame1(1)-frame2(1))^2)+((frame1(2)-frame2(2))^2)); % distance in cm
Dist_vect = Dist_vect / 100; % Distance in m
if isnan(Dist_vect)
    Distance(ii,1) = 0;
elseif ~isnan(Dist_vect)
    Distance(ii,1) = Dist_vect;
end

end

Timelapse = time(2) - time(1);

Speed = zeros(length(Distance),1);
for ii = 1:length(Distance)
Speed(ii,1) = Distance(ii,1)/Timelapse; % Speed in m/s
end
% plot(time(2:end),Speed)

%% 9_ Remove Jumpy speed

Jumpy = Speed > .8; % 0 not true, 1 true
JumpyFR = find(Jumpy == 1); % frames where Jumpy is 1

% plot track for the suspected jumpy detections
for ii = 1:length(JumpyFR)
A = track(JumpyFR(ii)-5:JumpyFR(ii)+5,1); 
B = track(JumpyFR(ii)-5:JumpyFR(ii)+5,2);
A2= track(JumpyFR(ii)-5:JumpyFR(ii)-1,1);
B2= track(JumpyFR(ii)-5:JumpyFR(ii)-1,2);
A3= track(JumpyFR(ii)+1:JumpyFR(ii)+5,1);
B3= track(JumpyFR(ii)+1:JumpyFR(ii)+5,2);

figure;
imshow(inArea)
hold on;
plot(A,B,'b','LineWidth',1);
plot(A2,B2,'b+');
plot(A3,B3,'g+');
plot(track(JumpyFR(ii),1),track(JumpyFR(ii),2),'ko');
legend('Track','PAST position','FUTURE position','PRESENT Position','Location','northeastoutside');
end

% ISjumpy = []; % You need it
ISjumpy = ones(length(JumpyFR),1); % Assumes that all is jumpy

%% 10_ Check the plots and make the ISjumpy variable
% ISjumpy is a 0 = not jumpy, 1 = jumpy.

for ii = 1:length(ISjumpy)
    if ISjumpy(ii) == 1
        C = JumpyFR(ii);
        Speed(C,1) = 0;
    elseif ISjumpy(ii) == 0
    end    
end

%% 11_ Calculations

plot(time(2:end),Speed);
int = time(numel(time))/(length(Speed)-1);
t2  = 0:int:time(numel(time));
% saveas(gcf, name(1:end-5), 'jpg'); % by Kinga
mSpeed = mean(Speed(t2>30)); % mean
mdSpeed = median(Speed(t2>30)); % median
sdSpeed = std(Speed(t2>30)); % standard dev
area = trapz(Speed(t2>30)); % area under curve

%% 12_ Save things

save([name(1:10) '.mat'], 'Speed', 'mSpeed', 'mdSpeed', 'sdSpeed', 'area', '-append');

%% 13_ Speed analysis
% Variables calculated and saved in a table: 
% - tone_Speed:     Speed during the tone;
% - pretone_Speed:  Speed before the tone starts; 
% - light_Speed:    Speed during the light; 
% - prelight_Speed: Speed before the light starts; 

% Seek the needed files
video_list = dir(['E','*.mp4']);
video_list = {video_list.name}';

% Set the name of the table
A = video_list{1};
totaltab = genvarname(A(5:end-26));
tab = table;

sessions = [3:4]; % Number of sessions to analyse
for kk = sessions 
    
    % Seek the needed files and load variables
    video = video_list{kk};
    Trackfile = [video(1:10),'.mat']; % Set file to save the time variable
    TTLfile   = [video(1:10),'_TTLs.mat']; % Set the file to save
    load(TTLfile,'TTL_light','TTL_sound','delay','TTL_sound_video');
    load(Trackfile,'Speed','times');
    delay = zeros(length(TTL_sound_video),1);
    
    try
        t2  = linspace(0,times(numel(times)),length(Speed));
    catch
        int = times(numel(times))/(length(Speed)-1);
        t2  = 0:int:times(numel(times));
    end
    
    % Generate tables
    SessionData = table; % Table for each session 
    tabname = video(1:end-26);

    % Synchronisation of TTLs to the video, delay calculated in step 2.
    sound = TTL_sound_video - delay;
%     sound = 150;
%     light = TTL_light - delay;

    for ii = 1:length(sound)
        SessionData.Track(ii,:) =[tabname '_s' num2str(ii,'%02.0f')];
            
        % Set the times of the TTL for each trial
        if exist('sound','var')
            sd = sound(ii)-0.2;
%             sk = sound(ii)-0.2+18;
        end
        
        if exist('light','var')
            li = light(ii)-0.2;
        end
        
        % Sound-related data ----------------------------------------------
        if exist('sd','var')
            % To extract speed the animal is moving at during the tone
            toneDur   = (t2>=sd & t2<=(sd+20)); toneDur(1,end) = 0; % hold on; plot((t2),toneDur*0.06,'g'); hold off
            toneSpeed = Speed(toneDur == 1);

            mToneSpeed  = mean(toneSpeed); % mean
            mdToneSpeed = median(toneSpeed); % median
            sdToneSpeed = std(toneSpeed); % standard dev
            areaToneSpeed = trapz(toneSpeed); % area under curve
            
            % To extract speed the animal is moving at 30s pretone
            pretoneDur   = (t2>=(sd-30) & t2<=sd); pretoneDur(1,end) = 0;
            pretoneSpeed = Speed(pretoneDur == 1);

            mPretoneSpeed  = mean(pretoneSpeed); % mean
            mdPretoneSpeed = median(pretoneSpeed); % median
            sdPretoneSpeed = std(pretoneSpeed); % standard dev
            areaPretoneSpeed = trapz(pretoneSpeed); % area under curve
            
            SessionData.tone_mSpeed(ii,:)    = mToneSpeed;
            SessionData.tone_mdSpeed(ii,:)   = mdToneSpeed;
            SessionData.tone_sdSpeed(ii,:)   = sdToneSpeed;
            SessionData.tone_areaSpeed(ii,:) = areaToneSpeed;
            SessionData.pretone_mSpeed(ii,:)    = mPretoneSpeed;
            SessionData.pretone_mdSpeed(ii,:)   = mdPretoneSpeed;
            SessionData.pretone_sdSpeed(ii,:)   = sdPretoneSpeed;
            SessionData.pretone_areaSpeed(ii,:) = areaPretoneSpeed;
        end

        % Light-related data ----------------------------------------------
        if exist('li','var')
            % To extract speed the animal is moving at during the light
            lightDur = (t2>=li & t2<=(li+10)); lightDur(1,end) = 0; % hold on; plot((t2),lightDur*0.06,'g'); hold off
            lightSpeed = Speed(lightDur == 1);

            mLightSpeed  = mean(lightSpeed); % mean
            mdLightSpeed = median(lightSpeed); % median
            sdLightSpeed = std(lightSpeed); % standard dev
            areaLightSpeed = trapz(lightSpeed); % area under curve

            % To extract speed the animal is moving at 30s prelight
            prelightDur = (t2>=(li-30) & t2<=li); prelightDur(1,end) = 0;
            prelightSpeed = Speed(lightDur == 1);

            mPrelightSpeed  = mean(prelightSpeed); % mean
            mdPrelightSpeed = median(prelightSpeed); % median
            sdPrelightSpeed = std(prelightSpeed); % standard dev
            areaPrelightSpeed = trapz(prelightSpeed); % area under curve
            
            SessionData.light_mSpeed(ii,:)    = mLightSpeed;
            SessionData.light_mdSpeed(ii,:)   = mdLightSpeed;
            SessionData.light_sdSpeed(ii,:)   = sdLightSpeed;
            SessionData.light_areaSpeed(ii,:) = areaLightSpeed;
            SessionData.prelight_mSpeed(ii,:)    = mPrelightSpeed;
            SessionData.prelight_mdSpeed(ii,:)   = mdPrelightSpeed;
            SessionData.prelight_sdSpeed(ii,:)   = sdPrelightSpeed;
            SessionData.prelight_areaSpeed(ii,:) = areaPrelightSpeed;
        end
        
        % Shock-related data ----------------------------------------------
        if exist('sk','var')
            % To extract speed the animal is moving at during the shock
            shockDur   = (t2>=sk & t2<=(sk+2)); shockDur(1,end) = 0; % hold on; plot((t2),toneDur*0.06,'g'); hold off
            shockSpeed = Speed(shockDur == 1);

            mShockSpeed  = mean(shockSpeed); % mean
            mdShockSpeed = median(shockSpeed); % median
            sdShockSpeed = std(shockSpeed); % standard dev
            areaShockSpeed = trapz(shockSpeed); % area under curve
            
            % To extract speed the animal is moving at 5s preshock
            preshockDur   = (t2>=(sk-5) & t2<=sk); preshockDur(1,end) = 0;
            preshockSpeed = Speed(preshockDur == 1);

            mPreshockSpeed  = mean(preshockSpeed); % mean
            mdPreshockSpeed = median(preshockSpeed); % median
            sdPreshockSpeed = std(preshockSpeed); % standard dev
            areaPreshockSpeed = trapz(preshockSpeed); % area under curve
            
            % To extract speed the animal is moving at 5s postshock
            postshockDur   = (t2>=(sk+2) & t2<=sk+7); postshockDur(1,end) = 0;
            postshockSpeed = Speed(postshockDur == 1);

            mPostshockSpeed  = mean(postshockSpeed); % mean
            mdPostshockSpeed = median(postshockSpeed); % median
            sdPostshockSpeed = std(postshockSpeed); % standard dev
            areaPostshockSpeed = trapz(postshockSpeed); % area under curve
            
            SessionData.shock_mSpeed(ii,:)    = mShockSpeed;
            SessionData.shock_mdSpeed(ii,:)   = mdShockSpeed;
            SessionData.shock_sdSpeed(ii,:)   = sdShockSpeed;
            SessionData.shock_areaSpeed(ii,:) = areaShockSpeed;
            SessionData.preshock_mSpeed(ii,:)    = mPreshockSpeed;
            SessionData.preshock_mdSpeed(ii,:)   = mdPreshockSpeed;
            SessionData.preshock_sdSpeed(ii,:)   = sdPreshockSpeed;
            SessionData.preshock_areaSpeed(ii,:) = areaPreshockSpeed;
            SessionData.postshock_mSpeed(ii,:)    = mPostshockSpeed;
            SessionData.postshock_mdSpeed(ii,:)   = mdPostshockSpeed;
            SessionData.postshock_sdSpeed(ii,:)   = sdPostshockSpeed;
            SessionData.postshock_areaSpeed(ii,:) = areaPostshockSpeed;
        end

    end

%      save([tabname '_analysis.mat'],'SessionData');
     tab = [tab; SessionData]; %#ok<AGROW>
     clearvars SessionData
end

eval([totaltab '= [tab];']);

%% If a second ROI is needed

% 2nd ROI
Plat = exist('inPlat','var'); % Idem as for 'inArea'
if Plat == 1    
elseif Plat ~= 1
    video = name;
    videoSource = vision.VideoFileReader(video); % Read the video
    frameSize = size(step(videoSource));
    videoPlayer = vision.VideoPlayer('Position',[100 100 [frameSize(2),...
    frameSize(1)]+30]); 

    [X2,Y2] = meshgrid(1:frameSize(2),1:frameSize(1)); 
    plat = step(videoSource);
    figure('name','Plat'); imshow(step(videoSource));
    [x2,y2] = ginput;
    inPlat = inpolygon(X2(:),Y2(:),x2,y2);
    inPlat = double(reshape(inPlat,[frameSize(1) frameSize(2)]));
    plat(~inPlat) = 0;
    hold on; imshow(plat); hold off;
    close Plat;
end

% INOUT of the 2nd ROI. 1 == IN, 0 == OUT.
lick_inout = zeros(1,length(track));
track2(:,1) = nan2zero(track(:,1));
track2(:,2) = nan2zero(track(:,2));
for ii = 1:length(track)
    position = zeros(frameSize(1),frameSize(2));
    position(round(track2(ii,2)+1),round(track2(ii,1))+1)= 1;
    INplat = repmat(position & inPlat,[1 1 1]);
    position(INplat == 0) = [];

    if isempty(position)==0
         lick_inout(ii) = 1; % The mice is IN the platform
    elseif isempty(position)==1
         lick_inout(ii) = 0; % The mice is OUT the platform
    end
end

%% 14_ Put together the animal data in the desired groups - Control
Ctrl = struct;

parentpath = pwd;
childpath = dir('Str*');
childpath = {childpath.name}';

% My inhibited (TelC) animals are #1, 3, 7, 8. Inhibited = [1, 3, 7, 8];
% My control (Control) animals are #2, 4, 5, 6. Control = [2, 4, 5, 6];
Control = [2, 4, 5, 6];
for ii = Control
    cd([parentpath,filesep,childpath{ii}])
    Animal = dir(['Str*', num2str(ii), '*Analysis.mat']);
    Animal = {Animal.name}';
    load(Animal{1,:});
    
    s = whos('*FCExt2min*'); 
    Tablename = s.name;
    Data = evalin('base',[Tablename '(:,:);']);
    
    VarName = Data.Properties.VariableNames;
    
    if ii == Control(1)
        Ctrl.(VarName{2}) = Data.(VarName{2});     
        Ctrl.(VarName{3}) = Data.(VarName{3}); 
        Ctrl.(VarName{5}) = Data.(VarName{5});     
        Ctrl.(VarName{6}) = Data.(VarName{6}); 
        Ctrl.(VarName{7}) = Data.(VarName{7});     
        Ctrl.(VarName{9}) = Data.(VarName{9}); 
%         Ctrl.(VarName{10}) = Data.(VarName{10}); 
%         Ctrl.(VarName{11}) = Data.(VarName{11});     
%         Ctrl.(VarName{13}) = Data.(VarName{13}); 
%         Ctrl.(VarName{14}) = Data.(VarName{14}); 
%         Ctrl.(VarName{15}) = Data.(VarName{15});     
%         Ctrl.(VarName{17}) = Data.(VarName{17});
%         Ctrl.(VarName{18}) = Data.(VarName{18});
%         Ctrl.(VarName{19}) = Data.(VarName{19});
%         Ctrl.(VarName{21}) = Data.(VarName{21});
    elseif ii ~= Control(1)
        Ctrl.(VarName{2}) = [Ctrl.(VarName{2}),Data.(VarName{2})];     
        Ctrl.(VarName{3}) = [Ctrl.(VarName{3}),Data.(VarName{3})]; 
        Ctrl.(VarName{5}) = [Ctrl.(VarName{5}),Data.(VarName{5})];     
        Ctrl.(VarName{6}) = [Ctrl.(VarName{6}),Data.(VarName{6})]; 
        Ctrl.(VarName{7}) = [Ctrl.(VarName{7}),Data.(VarName{7})];     
        Ctrl.(VarName{9}) = [Ctrl.(VarName{9}),Data.(VarName{9})]; 
%         Ctrl.(VarName{10}) = [Ctrl.(VarName{10}),Data.(VarName{10})]; 
%         Ctrl.(VarName{11}) = [Ctrl.(VarName{11}),Data.(VarName{11})];     
%         Ctrl.(VarName{13}) = [Ctrl.(VarName{13}),Data.(VarName{13})]; 
%         Ctrl.(VarName{14}) = [Ctrl.(VarName{14}),Data.(VarName{14})]; 
%         Ctrl.(VarName{15}) = [Ctrl.(VarName{15}),Data.(VarName{15})];     
%         Ctrl.(VarName{17}) = [Ctrl.(VarName{17}),Data.(VarName{17})]; 
%         Ctrl.(VarName{18}) = [Ctrl.(VarName{18}),Data.(VarName{18})]; 
%         Ctrl.(VarName{19}) = [Ctrl.(VarName{19}),Data.(VarName{19})]; 
%         Ctrl.(VarName{21}) = [Ctrl.(VarName{21}),Data.(VarName{21})]; 
    end
    
    clearvars -except Ctrl TelC parentpath childpath Control Inhibited
    cd(parentpath)
end

%% 15_ Put together the animal data in the desired groups - TelC
TelC = struct;

parentpath = pwd;
childpath = dir('Str*');
childpath = {childpath.name}';

% My inhibited (TelC) animals are #1, 3, 7, 8. Inhibited = [1, 3, 7, 8];
% My control (Control) animals are #2, 4, 5, 6. Control = [2, 4, 5, 6];
Inhibited = [1, 3, 7, 8];
for ii = Inhibited
    cd([parentpath,filesep,childpath{ii}])
    Animal = dir(['Str*', num2str(ii), '*Analysis.mat']);
    Animal = {Animal.name}';
    load(Animal{1,:});
    
    s = whos('*FCExt2min*'); 
    Tablename = s.name;
    Data = evalin('base',[Tablename '(:,:);']);
    
    VarName = Data.Properties.VariableNames;
    
    if ii == Inhibited(1)
        TelC.(VarName{2}) = Data.(VarName{2});     
        TelC.(VarName{3}) = Data.(VarName{3}); 
        TelC.(VarName{5}) = Data.(VarName{5});     
        TelC.(VarName{6}) = Data.(VarName{6}); 
        TelC.(VarName{7}) = Data.(VarName{7});     
        TelC.(VarName{9}) = Data.(VarName{9}); 
%         TelC.(VarName{10}) = Data.(VarName{10}); 
%         TelC.(VarName{11}) = Data.(VarName{11});     
%         TelC.(VarName{13}) = Data.(VarName{13}); 
%         TelC.(VarName{14}) = Data.(VarName{14}); 
%         TelC.(VarName{15}) = Data.(VarName{15});     
%         TelC.(VarName{17}) = Data.(VarName{17});
    elseif ii ~= Inhibited(1)
        TelC.(VarName{2}) = [TelC.(VarName{2}),Data.(VarName{2})];     
        TelC.(VarName{3}) = [TelC.(VarName{3}),Data.(VarName{3})]; 
        TelC.(VarName{5}) = [TelC.(VarName{5}),Data.(VarName{5})];     
        TelC.(VarName{6}) = [TelC.(VarName{6}),Data.(VarName{6})]; 
        TelC.(VarName{7}) = [TelC.(VarName{7}),Data.(VarName{7})];     
        TelC.(VarName{9}) = [TelC.(VarName{9}),Data.(VarName{9})]; 
%         TelC.(VarName{10}) = [TelC.(VarName{10}),Data.(VarName{10})]; 
%         TelC.(VarName{11}) = [TelC.(VarName{11}),Data.(VarName{11})];     
%         TelC.(VarName{13}) = [TelC.(VarName{13}),Data.(VarName{13})]; 
%         TelC.(VarName{14}) = [TelC.(VarName{14}),Data.(VarName{14})]; 
%         TelC.(VarName{15}) = [TelC.(VarName{15}),Data.(VarName{15})];     
%         TelC.(VarName{17}) = [TelC.(VarName{17}),Data.(VarName{17})]; 
    end
    
    clearvars -except Ctrl TelC parentpath childpath Control Inhibited
    cd(parentpath)
end

%% 16_ Statistical analysis

DataT = TelC3an_LickB.Lick_all;
DataC = Ctrl3an_LickB.Lick_all;
name = 'Analysis3an_LickB_all';
kk = 7; % Set to 3 for session, to 9 for day analysis

p_val = zeros(2,1);
C_dat = zeros(2,2);
T_dat = zeros(2,2);

for ii = 1:kk:size(DataT,1)
    
mT = nanmean(DataT(ii:ii+(kk-1),:),2);
mC = nanmean(DataC(ii:ii+(kk-1),:),2);
[h,p] = ttest2(nanmean(DataT(ii:ii+(kk-1),:),2),nanmean(DataC(ii:ii+(kk-1),:),2));
% % % For comparisson of the individual trials
% mT = nanmean(DataT(ii,:),2);
% mC = nanmean(DataC(ii,:),2);
% [h,p] = ttest2(nanmean(DataT(ii,:),2),nanmean(DataC(ii,:),2));

n = ceil(ii/kk);
p_val(n,1) = p;
T_dat(n,1) = nanmean(mT);
T_dat(n,2) = nanstd(mT);
C_dat(n,1) = nanmean(mC);
C_dat(n,2) = nanstd(mC);
% For comparisson of the individual trials
% T_dat(n,1) = nanmean(DataT(ii,:),2);
% T_dat(n,2) = nanstd(DataT(ii,:));
% C_dat(n,1) = nanmean(DataC(ii,:),2);
% C_dat(n,2) = nanstd(DataC(ii,:));

end

errorbar(T_dat(:,1),T_dat(:,2),'-o','LineWidth',2,'MarkerFaceColor','auto');
xlim([0 length(T_dat(:,1))+1]);
hold on;
errorbar(C_dat(:,1),C_dat(:,2),'-o','LineWidth',2,'MarkerFaceColor','auto');
y = get(gca,'ylim');
pval = double(p_val<=0.05)*(round(y(2)+abs(y(1)),3));
pval(pval == 0)= NaN;
plot(pval(:,1),'k*');
legend('TelC','Control','p-val < 0.05','Location', 'northwest');
ylim([round(y(1),3) (round(y(2)+abs(y(1)),2)+round(abs(y(1)),3))+0.002]);
% xticks([1 2])
xlabel('Session #');
ylabel('Speed (ms^-^1)');
set(gca,'fontname','arial');
hold off;
saveas(gcf,[name],'png')

save(name,'T_dat','C_dat','p_val');
