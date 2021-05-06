%% TRAQUING ANALYSIS
% Provides the steps to analyse the results from tracking the animals in an
% Active avoidance paradigm or similar behaviour. The animal behaviour was
% recorded with a webcam in an .mp4 video and the licking was recorded
% trough the analog input to an Intan Evaluation Board (Intan). The TTLs
% were split to syncrhonise the behaviour of the animal with the stimuli
% and the licking, then recorded with the Intan Evaluation Board.
%
% To track the movement of the animals use the Tracking_detection.mat script.
% Generates the track and inout variables to be used in this script.
% 
% STEP 1: Calculates and saves the TTLs. We use SHOCK, SOUND and LIGHT.
% STEP 2: Calculates the TTLs in the video to synchronise the video and the
% intan recording. A) using the sound; B) using the light.
% STEP 3: Calculates the time the animals spent in platform at given
% conditions.
%
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2020
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
%
% MATLAB toolboxes: - Computer Vision System 
%                   - Image Processing 
% -------------------------------------------------------------------------
 
%% STEP 1: Save TTLs
% The load_ttl script requires to manually enter the variable to use to do
% the ttl detection. Just tipe: shock, sound or light when asked.

Session = 'S';
% Seek .rhd files to extract the TTLs
TTL_list = dir([Session,'*.rhd']); 
TTL_list = {TTL_list.name}';

% Shock and digital signals for TTLs --------------------------------------
for ii = 1:size(TTL_list,1)
    name  = TTL_list{ii};
    read_Intan_RHD2000_file(name,cd)
    light = board_dig_in_data(1,:);
    sound = board_dig_in_data(2,:);
    shock = board_dig_in_data(3,:);

    newname = [name(:,1:end-17), 'TTLs.mat'];
    save(newname,'light','sound','shock','t_dig');
    
    disp('Now calculating: Shock TTL');
    TTL_shock = load_ttl(newname,'origin','IN','mPD',2,'mPW',0.005);
    save(newname,'TTL_shock','-append');
    
    clearvars light sound shock TTL_shock
end
disp('FINISHED! Shock TTL');

% Seek newly saved files
TTL_list = dir([Session,'*TTLs.mat']);
TTL_list = {TTL_list.name}';

% Sound -------------------------------------------------------------------
for ii = 1:size(TTL_list,1)
    TTLfile = TTL_list{ii};
    load(TTLfile,'sound');
    
    disp('Now calculating: Sound TTL');
    TTL_sound = load_ttl(TTLfile,'origin','IN','mPD',2,'mPW',0.01);
    save(TTLfile,'TTL_sound','-append');
    
    clearvars sound TTL_sound
end
disp('FINISHED! Sound TTL');

% Light -------------------------------------------------------------------
for ii = 1:size(TTL_list,1)
    TTLfile = TTL_list{ii};
    load(TTLfile,'light');
    
    disp('Now calculating: Light TTL');
    TTL_light = load_ttl(TTLfile,'origin','IN','mPD',2,'mPW',0.01);
    save(TTLfile,'TTL_light','-append');
    
    clearvars light TTL_light
end
disp('FINISHED! Light TTL');

%% STEP 2A: Synchronisation of intan TTLs with the video. 
% Uses the start of the sound in the video to set the time difference with
% the sound_TTL: delay.

freq = [10900 11100]; % 10990 11010 Frequencies for the bandpass. Depend on the sound used.
fig  = 1; % To plot the filtered signal and check that the detection was good. 
nSD  = 2.5; % Num times to multiply the SD for the detection. Check with plot.

% Seek the videos  
video_list = dir(['S' '*.mp4']); 
video_list = {video_list.name}';

for ii = 1:length(video_list)
    video = video_list{ii};
    Trackfile = [video(1:8),'_track.mat']; % Set file to save the time variable
    TTLfile   = [video(1:8),'_TTLs.mat']; % Set the file to save
    
    info  = audioinfo(video);
    audio = audioread(video);
    times = linspace(0,info.Duration,length(audio));
    fs    = 48000; % Audio sampling frequency
    wn    = [freq(1) freq(2)]/(fs/2); % Bandpass sound 
    [b,a] = butter(3,wn);
    f     = filter(b,a,audio(:,1));
    fm    = mean(f);
    fsd   = std(f);
    ft    = times > 105;
    fth   = f > fm+fsd*nSD;
    
    figure;
    findpeaks(double(fth(ft)),times(ft),'MinPeakDistance',30);
    [peaks,locs] = findpeaks(double(fth(ft)),times(ft),'MinPeakDistance',30); % ,'MinPeakWidth',5
    TTL_sound_video  = locs;
            
    if fig == 1
        hold on;
        plot(times, f*200)
        plot([0 max(times)],[fm+fsd*nSD*200 fm+fsd*nSD*200]);
        hold off;
    elseif fig == 0
        close figure
    end
    
    save(Trackfile,'times','-append');
    
    load(TTLfile,'TTL_sound');
    try 
        delay = (TTL_sound - TTL_sound_video);
        allD(ii,:) = delay;
        save(TTLfile,'delay','TTL_sound_video','-append');
    catch
        disp(['UUPS! Something went wrong in ',video(1:10),'. Please review!']);
        disp(['The number of TTLs is ', num2mstr(length(TTL_sound_video))]);
    end
    
end

%% STEP 2B: Synchronisation of intan TTLs with the video. 
% Uses the start of the light in the video to set the time difference with
% the light_TTL: delay.

% Use the tracking sript and set the ROI in the sourse of light, thus
% generating a tracking of the light that gives values when the light is
% on. Then use the Track variable, see below.

% Seek the needed files
video_list = dir(['S','*.mp4']);
video_list = {video_list.name}';

for ii = 1:length(video_list)
    video = video_list{ii};
    Trackfile = [video(1:10),'.mat']; % Set file to save the time variable
    TTLfile   = [video(1:10),'_TTLs.mat']; % Set the file to save
    
    load(Trackfile);
    info = audioinfo(video);
    times = linspace(0,info.Duration,length(light_inout));
    times = transpose(times);
    save(Trackfile,'times','-append');

    TTL_light_video = load_ttl(Trackfile,'origin','WS','n',60,'mPP',0.5,'mPD',0.1,'mPW',9,'ttl_dur',11);
    TTL_light_video = transpose(TTL_light_video);
    
    load(TTLfile,'TTL_light');
    try
        delay = (TTL_light - TTL_light_video);
        disp(['Delay of ', video(1:10), ' is ', num2str(delay)]);
        save(TTLfile,'delay','TTL_light_video','-append');
    catch
        disp(['UUPS! Something went wrong in ',video(1:10),'. Please review!']);
        disp(['The number of TTLs is ', num2mstr(length(TTL_light_video))]);
    end
    
%     figure;
%     plot(times,light_inout);
%     set(gcf,'Position',[5,672,1914,306]);
%     hold on;
%     plot(TTL_light_video,ones(1,length(TTL_light_video)),'x');
%     title(['Session ' video(1:10)]);
%     plot(TTL_light,ones(1,length(TTL_light)),'+');
%     hold off;
%     clearvars TTL_light_video TTL_light times light_inout delay
    
end

%% STEP 3: Calculate the time spent in platform (%).
% Variables calculated and saved in a table: 
% - tIN_tone:     Time IN the platform during the tone;
% - tIN_pretone:  Time IN the platform before the tone starts; 
% - tIN_light:    Time IN the platform during the light; 
% - tIN_prelight: Time IN the platform before the light starts; 
% - tIN_all:      Time IN the platform; 
% - tIN_shock:    Shock result; 
% - tonemin_lat:  First time to enter the platform after the start of the
%                 tone; 
% - tonemax_lat:  Last time to enter the platform after the start of the
%                 tone; 
% - lightmin_lat: First time to enter the platform after the start of the
%                 light;
% - lightmax_lat: Last time to enter the platform after the start of the
%                 light.

% Seek the needed files
video_list = dir(['*PrL2','*.mp4']);
video_list = {video_list.name}';

% To plot the figures
fig = 1;

% Set the name of the table
A = video_list{1};
totaltab = genvarname(A(5:end-26));
tab = table;

sessions = [1:10]; % Number of sessions to analyse
for kk = sessions 
    
    % Seek the needed files and load variables
    video = video_list{kk};
    Trackfile = [video(1:8),'_track.mat']; % Set file to save the time variable
    TTLfile   = [video(1:8),'_TTLs.mat']; % Set the file to save
%     if kk == 11 || kk == 14
%         Trackfile = [video(1:8),'_Test_Track.mat']; % Set file to save the time variable
%         TTLfile   = [video(1:8),'_Test_TTLs.mat']; % Set the file to save
%     end
    load(TTLfile,'TTL_shock','TTL_light','TTL_sound','delay');
    load(Trackfile,'inout','times');
    
    try
        t2  = linspace(0,times(numel(times)),length(inout));
    catch
        int = times(numel(times))/(length(inout)-1);
        t2  = 0:int:times(numel(times));
    end

    % Generate tables
    SessionData = table; % Table for each session 
    tabname = video(1:end-26);
    SessionData.Track =[[tabname '_s1'];[tabname '_s2'];[tabname '_s3'];...
        [tabname '_s4'];[tabname '_s5'];[tabname '_s6'];...
        [tabname '_s7']];

    % Synchronisation of TTLs to the video, delay calculated in step 2.
    shock = TTL_shock - delay; 
    sound = TTL_sound - delay;
    light = TTL_light - delay;

    for ii = 1:length(shock)
        
        % Set the times of the TTL for each trial
        sk = shock(ii)-0.2; 
        sd = sound(ii)-0.2;
        li = light(ii)-0.2;
%         sk = floor(round(shock(ii)*10,1))/10; % Shock time up to 1st decimal
%         sd = floor(round(sound(ii)*10,1))/10;
%         li = floor(round(light(ii)*10,1))/10;

        % Remove NaN and fill with IN values (1) or OUT values (0)
        idx = find(isnan(inout));
        inout2 = inout;
        if isempty(idx) == 1
        elseif idx(1,1) == 1
            inout2(1,1) = 0;
            idx(1,:) = []; 
            for i=1:length(idx)
            inout2(idx) = inout2(idx-1);
            end
        elseif idx(1,1) ~= 1
            for i=1:length(idx)
            inout2(idx) = inout2(idx-1);
            end
        end
                
        % To extract time the animal is inout the platform during the tone
        toneDur = (t2>=sd & t2<=(sd+21)); toneDur(1,end) = 0; % hold on; plot((t2),toneDur*0.06,'g'); hold off
        toneIO = transpose(inout2) & toneDur;
        toneIO(1,1) = 0;

        A = ischange(double(toneIO),'MaxNumChanges',length(toneIO));
        B = t2(A); 
        if mod(length(B),2)==0 
        elseif mod(length(B),2)~=0
            B = [B,t2(end)];
        end
        tone_t = transpose(B(:,1:2:end));  % odd matrix
        tone_t(:,2) = transpose(B(:,2:2:end)); % even matrix

        A = ischange(double(toneDur),'MaxNumChanges',length(toneDur));
        tone_d = t2(A);

        % And the latency to enter the platform during tone
        if isempty(B)
            tonemin_lat = 21;
            tonemax_lat = 21;
        elseif ~isempty(B)
            tonemin_lat = tone_t(1,1) - tone_d(1,1);
            tonemax_lat = tone_t(end,1) - tone_d(1,1);
        end

        % To extract time the animal is inout the platform during the light
        lightDur = (t2>=li & t2<=(li+11)); lightDur(1,end) = 0; % hold on; plot((t2),lightDur*0.06,'g'); hold off
        lightIO = transpose(inout2) & lightDur;
        lightIO(1,1) = 0;

        A = ischange(double(lightIO),'MaxNumChanges',length(lightIO));
        B = t2(A); 
        if mod(length(B),2)==0 
        elseif mod(length(B),2)~=0
            B = [B,t2(end)];
        end
        light_t = transpose(B(:,1:2:end));  % odd matrix
        light_t(:,2) = transpose(B(:,2:2:end)); % even matrix

        A = ischange(double(lightDur),'MaxNumChanges',length(lightDur));
        light_d = t2(A);

        % And the latency to enter the platform during light
        if isempty(B)
            lightmin_lat = 11;
            lightmax_lat = 11;
        elseif ~isempty(B)
            lightmin_lat = light_t(1,1) - light_d(1,1);
            lightmax_lat = light_t(end,1) - light_d(1,1);
        end

        % To extract time the animal is inout the platform 2min pretone
        pretoneDur = (t2>=(sd-90) & t2<=sd); pretoneDur(1,end) = 0;
        pretoneIO = transpose(inout2) & pretoneDur;
        pretoneIO(1,1) = 0;

        A = ischange(double(pretoneIO),'MaxNumChanges',length(pretoneIO));
        B = t2(A); 
        if mod(length(B),2)==0 
        elseif mod(length(B),2)~=0
            B = [B,t2(end)];
        end
        pretone_t = transpose(B(:,1:2:end));  % odd matrix
        pretone_t(:,2) = transpose(B(:,2:2:end)); % even matrix

        A = ischange(double(pretoneDur),'MaxNumChanges',length(pretoneDur));
        pretone_d = t2(A);

        % To extract time the animal is inout the platform 2min prelight
        prelightDur = (t2>=(li-90) & t2<=li); prelightDur(1,end) = 0;
        prelightIO = transpose(inout2) & prelightDur;
        prelightIO(1,1) = 0;

        A = ischange(double(prelightIO),'MaxNumChanges',length(prelightIO));
        B = t2(A); 
        if mod(length(B),2)==0 
        elseif mod(length(B),2)~=0
            B = [B,t2(end)];
        end
        prelight_t = transpose(B(:,1:2:end));  % odd matrix
        prelight_t(:,2) = transpose(B(:,2:2:end)); % even matrix

        A = ischange(double(prelightDur),'MaxNumChanges',length(prelightDur));
        prelight_d = t2(A);

        % To extract time the animal is inout the platform all
        allDur = (t2>=(sd-90) & t2<=(sd+21)); allDur(1,end) = 0;
        allIO = transpose(inout2) & allDur;
        allIO(1,1) = 0;

        A = ischange(double(allIO),'MaxNumChanges',length(allIO));
        B = t2(A); 
        if mod(length(B),2)==0 
        elseif mod(length(B),2)~=0
            B = [B,t2(end)];
        end
        all_t = transpose(B(:,1:2:end));  % odd matrix
        all_t(:,2) = transpose(B(:,2:2:end)); % even matrix

        A = ischange(double(allDur),'MaxNumChanges',length(allDur));
        all_d = t2(A);

        % To know if the animal gets shocked
        shockDur = (t2>=sk & t2<=(sk+3)); % hold on; plot((t2),shockDur*0.06,'g'); hold off
        shockIO = transpose(inout2) & shockDur;
        S = find(shockDur == 1);
        if (sum(abs(shockIO(S)-1))) >= 1
            shocked = 1;
        else 
            shocked = 0;
        end 

        % Calculations
        if isempty(tone_t)== 1
            tone_t = [0 0];
        end
        if isempty(tone_d) == 1
            tone_d = [0 0];
        end
        if isempty(pretone_t)== 1
            pretone_t = [0 0];
        end
        if isempty(light_t)== 1
            light_t = [0 0];
        end
        if isempty(prelight_t)== 1
            prelight_t = [0 0];
        end
        if isempty(all_t)== 1
            all_t = [0 0];
        end 

        SessionData.tIN_tone(ii,:) = 100*sum(tone_t(:,2)-tone_t(:,1))/(tone_d(:,2)...
            -tone_d(:,1));
        SessionData.tIN_pretone(ii,:) = 100*sum(pretone_t(:,2)-pretone_t(:,1))...
            /(pretone_d(:,2)-pretone_d(:,1));
        SessionData.tIN_light(ii,:) = 100*sum(light_t(:,2)-light_t(:,1))...
            /(light_d(:,2)-light_d(:,1));
        SessionData.tIN_prelight(ii,:) = 100*sum(prelight_t(:,2)-prelight_t(:,1))...
            /(prelight_d(:,2)-prelight_d(:,1));
        SessionData.tIN_all(ii,:) = 100*sum(all_t(:,2)-all_t(:,1))/(all_d(:,2)-all_d(:,1));
        SessionData.tIN_shock(ii,:) = shocked;
        SessionData.tonemin_lat(ii,:) = tonemin_lat;
        SessionData.tonemax_lat(ii,:) = tonemax_lat;
        SessionData.lightmin_lat(ii,:) = lightmin_lat;
        SessionData.lightmax_lat(ii,:) = lightmax_lat;
        
        if fig == 1
            area(([sk sk+3]),([1.3 1.3]),'FaceAlpha',0.5,'FaceColor',...
                [0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840]);
            hold on;
            area(([li li+11]),([1.2 1.2]),'FaceAlpha',0.5,'FaceColor',...
                [0.9290 0.6940 0.1250],'EdgeColor',[0.9290 0.6940 0.1250]);
            area(([sd sd+21]),([1.1 1.1]),'FaceAlpha',0.5,'FaceColor',...
                [0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410]);
            xlim([sd-2 sd+22]); 
            plot(t2(allDur), allIO(allDur),'k','LineWidth',2);
            ylim([0 1.5]);
            legend('Shock','Light','Sound','Time in platform','Location','northwest');
            hold off;
            saveas(gcf, [tabname '_' num2str(ii)], 'png');
            close gcf
        end

    end

%      save([tabname '_analysis.mat'],'SessionData');
     tab = [tab; SessionData];
     clearvars SessionData TTL_light TTL_shock TTL_sound delay
end

eval([totaltab '= [tab];']);

%% Alternative ploting of the track
for nn = 1:10
    Animal = dir(['S*', num2str(nn,'%02.0f'),'*PrL2', '*.png']);
    Animal = {Animal.name}';

    for ii = 1:length(Animal)
        A = imread(Animal{ii});
        h = subplot(2,4,ii);
        pos = get(h, 'Position');
        if ii <= 4
            pos(2) = 0.52;
            pos(3) = 0.19;
            pos(4) = 0.4;
            set(h, 'Position', pos );
        elseif ii > 4
            pos(2) = 0.05;
            pos(3) = 0.19;
            pos(4) = 0.4;
            set(h, 'Position', pos );
        end
        image(A)
        set(gca,'XColor', 'none','YColor','none')
        title(['Trial #', num2str(ii)]);
    end
    set(gcf,'Position',[205,293,1715,650]);
    name = Animal{1};
    saveas(gcf,name(1:8),'png')
    close all
end

%% STEP 4: Ploting the data

% newtab = load('BA2b_table.mat');
% newtab = struct2cell(newtab);
% Data = newtab{1,1};
Data = PrL2_lick;
Animal_Data = struct;
nn = 1; % Grouping number, to analyse data in groups of nn trials

Animal_Data.tone     = zeros(length(Data.Track)/nn,1);
Animal_Data.pretone  = zeros(length(Data.Track)/nn,1);
Animal_Data.light    = zeros(length(Data.Track)/nn,1);
Animal_Data.prelight = zeros(length(Data.Track)/nn,1);
Animal_Data.all      = zeros(length(Data.Track)/nn,1);
Animal_Data.shock    = zeros(length(Data.Track)/nn,1);
Animal_Data.tonemin_lat  = zeros(length(Data.Track)/nn,1);
Animal_Data.lightmin_lat = zeros(length(Data.Track)/nn,1);
Animal_Data.tonemax_lat  = zeros(length(Data.Track)/nn,1);
Animal_Data.lightmax_lat = zeros(length(Data.Track)/nn,1);
for ii = 1:nn:length(Data.Track)
    tone     = nanmean(Data.tIN_tone(ii:ii+(nn-1),:));
    pretone  = nanmean(Data.tIN_pretone(ii:ii+(nn-1),:));
    light    = nanmean(Data.tIN_light(ii:ii+(nn-1),:));
    prelight = nanmean(Data.tIN_prelight(ii:ii+(nn-1),:));
    all      = nanmean(Data.tIN_all(ii:ii+(nn-1),:));
    shock    = max(Data.tIN_shock(ii:ii+(nn-1),:));
    shock100 = nanmean(Data.tIN_shock(ii:ii+(nn-1),:))/100;
    tonemin_lat  = nanmean(Data.tonemin_lat(ii:ii+(nn-1),:));
    lightmin_lat = nanmean(Data.lightmin_lat(ii:ii+(nn-1),:));
    tonemax_lat  = nanmean(Data.tonemax_lat(ii:ii+(nn-1),:));
    lightmax_lat = nanmean(Data.lightmax_lat(ii:ii+(nn-1),:));
    
    n = ceil(ii/nn);
    Animal_Data.tone(n,1)     = tone;
    Animal_Data.pretone(n,1)  = pretone;
    Animal_Data.light(n,1)    = light;
    Animal_Data.prelight(n,1) = prelight;
    Animal_Data.all(n,1)      = all;
    Animal_Data.shock(n,1)    = shock;
    Animal_Data.shock100(n,1) = shock100;
    Animal_Data.tonemin_lat(n,1)  = tonemin_lat;
    Animal_Data.lightmin_lat(n,1) = lightmin_lat;
    Animal_Data.tonemax_lat(n,1)  = tonemax_lat;
    Animal_Data.lightmax_lat(n,1) = lightmax_lat;

end

figname = Data.Track(1,5:end-3);

subplot(3,1,1)
plot(Animal_Data.tone,'b-o','LineWidth',2,'MarkerFaceColor','b');
hold on;
plot(Animal_Data.pretone,'b--o','LineWidth',2,'MarkerFaceColor','b');
plot(Animal_Data.light,'c-o','LineWidth',2,'MarkerFaceColor','c');
plot(Animal_Data.prelight,'c--o','LineWidth',2,'MarkerFaceColor','c');
plot(Animal_Data.all,'r-o','LineWidth',2,'MarkerFaceColor','r');
hold off;
ylim([-5 105]);
xlim([0 length(Animal_Data.tone)+1]);
shock = Animal_Data.shock;
shock (shock == 0)= NaN;
hold on;
plot(shock,'o','MarkerEdgeColor',[1 0.41 0.16],'Markersize',10,'MarkerFaceColor',[1 0.41 0.16]);
legend('Time INplat - tone','Time INplat - pretone','Time INplat - light',...
    'Time INplat - prelight','Time INplat - total','Shock','Location', 'northeastoutside');
title('Time spent in the platform & Shocks received');
xlabel('Session #'); 
ylabel('% Of Time'); 
hold off;
% set(gca,'fontname','arial');  % Set font to arial
% saveas(gca,[figname, '_TimeIN'],'png');

subplot(3,1,2)
plot(Animal_Data.tonemin_lat,'k-o','LineWidth',2,'MarkerFaceColor','k');
hold on; plot(Animal_Data.tonemax_lat,'b-o','LineWidth',2,'MarkerFaceColor','b'); 
xlim([0 length(Animal_Data.tonemax_lat)+1]);
legend('Min latency (tone)','Max latency (tone)','Location', 'northeastoutside');
title('Latencies to enter the platform during the TONE');
xlabel('Session #'); 
ylabel('Time (s)');
hold off;
% set(gca,'fontname','arial');
% saveas(gca,[figname, '_ToneLat'],'png');

subplot(3,1,3)
plot(Animal_Data.lightmin_lat,'k-o','LineWidth',2,'MarkerFaceColor','k');
hold on; plot(Animal_Data.lightmax_lat,'c-o','LineWidth',2,'MarkerFaceColor','c'); 
xlim([0 length(Animal_Data.lightmax_lat)+1]);
legend('Min latency (light)','Max latency (light)','Location', 'northeastoutside');
title('Latencies to enter the platform during the LIGHT');
xlabel('Session #');
ylabel('Time (s)'); 
set(gca,'fontname','arial');
set(gcf,'Position',[471,154,1184,809]);
hold off;
% saveas(gca,[figname, '_LightLat'],'png');

% save([figname, '_Data'],'Animal_Data','-append');

%% STEP 5A: Put together the animal data in the desired groups - Inhibited

TelC = struct;
TelC.tone     = [];
TelC.pretone  = [];
TelC.light    = [];
TelC.prelight = [];
TelC.all      = [];
TelC.shock    = [];
TelC.tonemin_lat  = [];
TelC.lightmin_lat = [];
TelC.tonemax_lat  = [];
TelC.lightmax_lat = [];

parentpath = pwd;
childpath = dir('OF*');
childpath = {childpath.name}';

% My inhibited (TelC) animals are #1, 3, 7, 8
% My control (Control) animals are #2, 4, 5, 6
for ii = [5, 6]
    cd([parentpath,filesep,childpath{ii}])
    Animal = dir(['OF*', num2str(ii), '*_data.mat']);
    Animal = {Animal.name}';
    load(Animal{1,:});
    
    s = whos('OF*'); 
    Tablename = s.name;
    Data = evalin('base',[Tablename '(:,:);']);
    
    TelC.tone     = [TelC.tone, Data.tIN_tone];
    TelC.pretone  = [TelC.pretone, Data.tIN_pretone];
    TelC.light    = [TelC.light, Data.tIN_light];
    TelC.prelight = [TelC.prelight, Data.tIN_prelight];
    TelC.all      = [TelC.all, Data.tIN_all];
    TelC.shock    = [TelC.shock, Data.tIN_shock];
    TelC.tonemin_lat  = [TelC.tonemin_lat, Data.tonemin_lat];
    TelC.lightmin_lat = [TelC.lightmin_lat, Data.lightmin_lat];
    TelC.tonemax_lat  = [TelC.tonemax_lat, Data.tonemax_lat];
    TelC.lightmax_lat = [TelC.lightmax_lat, Data.lightmax_lat];
    
    clearvars -except Ctrl TelC parentpath childpath 
    cd(parentpath)
end

%% STEP 5B: Put together the animal data in the desired groups - Control

Ctrl = struct;
Ctrl.tone     = [];
Ctrl.pretone  = [];
Ctrl.light    = [];
Ctrl.prelight = [];
Ctrl.all      = [];
Ctrl.shock    = [];
Ctrl.tonemin_lat  = [];
Ctrl.lightmin_lat = [];
Ctrl.tonemax_lat  = [];
Ctrl.lightmax_lat = [];

parentpath = pwd;
childpath = dir('OF*');
childpath = {childpath.name}';

% My inhibited (TelC) animals are #1, 3, 7, 8
% My control (Control) animals are #2, 4, 5, 6
for ii = [1, 2, 4, 7]
    cd([parentpath,filesep,childpath{ii}])
    Animal = dir(['OF*', num2str(ii), '*_data.mat']);
    Animal = {Animal.name}';
    load(Animal{1,:});
    
    s = whos('OF*'); 
    Tablename = s.name;
    Data = evalin('base',[Tablename '(:,:);']);
    
    Ctrl.tone     = [Ctrl.tone, Data.tIN_tone];
    Ctrl.pretone  = [Ctrl.pretone, Data.tIN_pretone];
    Ctrl.light    = [Ctrl.light, Data.tIN_light];
    Ctrl.prelight = [Ctrl.prelight, Data.tIN_prelight];
    Ctrl.all      = [Ctrl.all, Data.tIN_all];
    Ctrl.shock    = [Ctrl.shock, Data.tIN_shock];
    Ctrl.tonemin_lat  = [Ctrl.tonemin_lat, Data.tonemin_lat];
    Ctrl.lightmin_lat = [Ctrl.lightmin_lat, Data.lightmin_lat];
    Ctrl.tonemax_lat  = [Ctrl.tonemax_lat, Data.tonemax_lat];
    Ctrl.lightmax_lat = [Ctrl.lightmax_lat, Data.lightmax_lat];
    
    clearvars -except Ctrl TelC parentpath childpath 
    cd(parentpath)
end

%% STEP 6: Statistical analysis

DataT = PrL1_lick.lick_all*100;
DataC = PrL2_lick.lick_all*100;
name = 'PrL_lick_all';
kk = 7; % Set to 3 for session, to 9 for day analysis

p_val = zeros(5,1);
C_dat = zeros(5,2);
T_dat = zeros(5,2);

for ii = 1:kk:length(DataC)
    
mT = nanmean(DataT(ii:ii+(kk-1),:),2);
mC = nanmean(DataC(ii:ii+(kk-1),:),2);
[h,p] = ttest2(nanmean(DataT(ii:ii+(kk-1),:),2),nanmean(DataC(ii:ii+(kk-1),:),2));

n = ceil(ii/kk);
p_val(n,1) = p;
T_dat(n,1) = nanmean(mT);
T_dat(n,2) = nanstd(mT);
C_dat(n,1) = nanmean(mC);
C_dat(n,2) = nanstd(mC);

end

pval = double(p_val<=0.07)*120;
pval(pval == 0)= NaN;
errorbar(T_dat(:,1),T_dat(:,2),'-o','LineWidth',2,'MarkerFaceColor','auto');
xlim([0 length(T_dat(:,1))+1]);
hold on;
errorbar(C_dat(:,1),C_dat(:,2),'-o','LineWidth',2,'MarkerFaceColor','auto');
plot(pval,'k*');
pval = double(p_val<=0.05)*120;
pval(pval == 0)= NaN;
plot(pval(:,1),'k*');
legend('PrL1','PrL2','p-val < 0.07','p-val < 0.05','Location', 'southeast');
ylim([-15 130]);
yticks([0 50 100])
xlabel('Session #');
ylabel('% Of Time');
set(gca,'fontname','arial');
hold off;
saveas(gcf,[name],'png')

save(name,'T_dat','C_dat','p_val');

%% Time spent on platform in 3 s bins during tone, Diehl et al 2018 eLife

video_list = dir(['*PrL2','*.mp4']);
video_list = {video_list.name}';

% Set the name of the table
A = video_list{1};
totaltab = [genvarname(A(5:end-26)),'_3sbins'];
tab = table;

sessions = [1:10]; % Number of sessions to analyse
for kk = sessions 
    % Seek the needed files and load variables
    video = video_list{kk};
    Trackfile = [video(1:8),'_track.mat']; % Set file to save the time variable
    TTLfile   = [video(1:8),'_TTLs.mat']; % Set the file to save
%     if kk == 11 || kk == 14
%         Trackfile = [video(1:7),'_Test_Track.mat']; % Set file to save the time variable
%         TTLfile   = [video(1:7),'_Test_TTLs.mat']; % Set the file to save
%     end
    load(TTLfile,'TTL_sound','delay');
    load(Trackfile,'inout','times');
    
    SessionData = table; % Table for each session 
    tabname = video(1:end-26);
    SessionData.Track =[[tabname '_s1'];[tabname '_s2'];[tabname '_s3'];...
        [tabname '_s4'];[tabname '_s5'];[tabname '_s6'];...
        [tabname '_s7']];
        
    try
        t2  = linspace(0,times(numel(times)),length(inout));
    catch
        int = times(numel(times))/(length(inout)-1);
        t2  = 0:int:times(numel(times));
    end
    
    % Synchronisation of TTLs to the video, delay calculated in step 2.
    sound = TTL_sound - delay;
    
    for ii = 1:length(sound)
        
        % Set the times of the TTL for each trial
        sd = sound(ii)-1;
        
        % Remove NaN and fill with IN values (1) or OUT values (0)
        idx = find(isnan(inout));
        inout2 = inout;
        if isempty(idx) == 1
        elseif idx(1,1) == 1
            inout2(1,1) = 0;
            idx(1,:) = []; 
            for i=1:length(idx)
            inout2(idx) = inout2(idx-1);
            end
        elseif idx(1,1) ~= 1
            for i=1:length(idx)
            inout2(idx) = inout2(idx-1);
            end
        end
                
        % To extract time the animal is inout the platform during the tone
        % 3min bins
        timeBin = t2>=(sd-3) & t2<=sd;
        toneIO = transpose(inout2) & timeBin;
        toneIO(1,1) = 0;

        A = ischange(double(toneIO),'MaxNumChanges',length(toneIO));
        B = t2(A); 
        if mod(length(B),2)==0 
        elseif mod(length(B),2)~=0
            B = [B,t2(end)];
        end
        tone_t = transpose(B(:,1:2:end));  % odd matrix
        tone_t(:,2) = transpose(B(:,2:2:end)); % even matrix

        A = ischange(double(timeBin),'MaxNumChanges',length(timeBin));
        tone_d = t2(A);

        if isempty(tone_t)== 1
            tone_t = [0 0];
        end
        if length(tone_d) == 1
            timeBin(end) = 0;
            A = ischange(double(timeBin),'MaxNumChanges',length(timeBin));
            tone_d = t2(A);
        elseif length(tone_d) < 1
                    tone_d = [0 0];
        end 

        SessionData.Bin0(ii) = 100*sum(tone_t(:,2)-tone_t(:,1))/(tone_d(:,2)...
            -tone_d(:,1));

        for jj = 1:7
            if     jj < 7
                timeBin = t2>=(sd+(jj-1)*3) & t2<=(sd+(jj-1)*3+3);
                toneIO = transpose(inout2) & timeBin;
                toneIO(1,1) = 0;

                A = ischange(double(toneIO),'MaxNumChanges',length(toneIO));
                B = t2(A); 
                if mod(length(B),2)==0 
                elseif mod(length(B),2)~=0
                    B = [B,t2(end)];
                end
                tone_t = transpose(B(:,1:2:end));  % odd matrix
                tone_t(:,2) = transpose(B(:,2:2:end)); % even matrix
                
                A = ischange(double(timeBin),'MaxNumChanges',length(timeBin));
                tone_d = t2(A);
        
                if isempty(tone_t)== 1
                    tone_t = [0 0];
                end
                if length(tone_d) == 1
                    timeBin(end) = 0;
                    A = ischange(double(timeBin),'MaxNumChanges',length(timeBin));
                    tone_d = t2(A);
                elseif length(tone_d) < 1
                    tone_d = [0 0];
                end 
                
                SessionData.(['Bin',num2str(jj)])(ii) = 100*sum(tone_t(:,2)-tone_t(:,1))/(tone_d(:,2)...
                    -tone_d(:,1));
            elseif jj == 7
                timeBin = t2>=(sd+(jj-1)*3) & t2<=(sd+(jj-1)*3+2);
                toneIO = transpose(inout2) & timeBin;
                toneIO(1,1) = 0;

                A = ischange(double(toneIO),'MaxNumChanges',length(toneIO));
                B = t2(A); 
                if mod(length(B),2)==0 
                elseif mod(length(B),2)~=0
                    B = [B,t2(end)];
                end
                tone_t = transpose(B(:,1:2:end));  % odd matrix
                tone_t(:,2) = transpose(B(:,2:2:end)); % even matrix
                
                A = ischange(double(timeBin),'MaxNumChanges',length(timeBin));
                tone_d = t2(A);
                
                if isempty(tone_t)== 1
                    tone_t = [0 0];
                end
                if length(tone_d) == 1
                    timeBin(end) = 0;
                    A = ischange(double(timeBin),'MaxNumChanges',length(timeBin));
                    tone_d = t2(A);
                elseif length(tone_d) < 1
                    tone_d = [0 0];
                end 
                
                SessionData.(['Bin',num2str(jj)])(ii) = 100*sum(tone_t(:,2)-tone_t(:,1))/(tone_d(:,2)...
                    -tone_d(:,1));
            end
        end
    end
    
    tab = [tab; SessionData];
    clearvars SessionData
end

eval([totaltab '= [tab];']);
% save([tabname '_analysis.mat'],'SessionData');

%% Statistical analysis 3bins
% # Session: 1 = (1:7), 2 = (8:15)...
S = 10;
b = S*7;
a = S*7 - 6;
DataT = PrL1_3sbins;
DataC = PrL2_3sbins;

Data_TelC(1:7,1) =  nanmean(DataT.Bin0(a:b,:),2);
Data_TelC(1:7,2) =  nanmean(DataT.Bin1(a:b,:),2);
Data_TelC(1:7,3) =  nanmean(DataT.Bin2(a:b,:),2);
Data_TelC(1:7,4) =  nanmean(DataT.Bin3(a:b,:),2);
Data_TelC(1:7,5) =  nanmean(DataT.Bin4(a:b,:),2);
Data_TelC(1:7,6) =  nanmean(DataT.Bin5(a:b,:),2);
Data_TelC(1:7,7) =  nanmean(DataT.Bin6(a:b,:),2);
Data_TelC(1:7,8) =  nanmean(DataT.Bin7(a:b,:),2);
Data_Ctrl(1:7,1) =  nanmean(DataC.Bin0(a:b,:),2);
Data_Ctrl(1:7,2) =  nanmean(DataC.Bin1(a:b,:),2);
Data_Ctrl(1:7,3) =  nanmean(DataC.Bin2(a:b,:),2);
Data_Ctrl(1:7,4) =  nanmean(DataC.Bin3(a:b,:),2);
Data_Ctrl(1:7,5) =  nanmean(DataC.Bin4(a:b,:),2);
Data_Ctrl(1:7,6) =  nanmean(DataC.Bin5(a:b,:),2);
Data_Ctrl(1:7,7) =  nanmean(DataC.Bin6(a:b,:),2);
Data_Ctrl(1:7,8) =  nanmean(DataC.Bin7(a:b,:),2);

for ii = 1:8
    [h,p] = ttest2(Data_TelC(:,ii),Data_Ctrl(:,ii));
    p_val(ii) = p;
end

pval = double(p_val<=0.07)*120;
pval(pval == 0)= NaN;
errorbar(nanmean(Data_TelC,1),nanstd(Data_TelC,1),'-o','LineWidth',2,'MarkerFaceColor','auto');
hold on;
errorbar(nanmean(Data_Ctrl,1),nanstd(Data_Ctrl,1),'-o','LineWidth',2,'MarkerFaceColor','auto');
plot(pval,'r*');
pval = double(p_val<=0.05)*120;
pval(pval == 0)= NaN;
plot(pval,'k*');
ylabel('Time IN platform (%)')
xlabel('Time (s)')
xticklabels([-3 0 3 6 9 12 15 18]);
ylim([0 130]);
xlim([0 9]);
title(['Mean time IN platform - S',num2str(S,'%02.0f')]);
lgd = legend('PrL1','PrL2','pval < 0.07','pval < 0.05','Location', 'southeast');
lgd.AutoUpdate = 'off';
plot([2 8],[125 125],'Color',[0.6 0.6 0.6],'LineWidth',10);
text(4.5,125.2,'TONE','FontWeight','bold','Color','k');
set(gca,'fontname','arial');
hold off;

saveas(gcf,['TimeINPlat_3sbins_S',num2str(S,'%02.0f')],'png')

%% Other plots
plot(Control_3.tone,'-o','LineWidth',2,'MarkerFaceColor','auto');
hold on;
shock = Control_3.shock;
shock (shock == 0)= NaN;
hold on;
plot(shock,'o','MarkerEdgeColor',[1 0.41 0.16],'Markersize',10,'MarkerFaceColor','y');
legend('Time in platform','Shock');
hold off;
ylim([-5 105]);
xlim([0 length(Control_3.tone)+1]);

%% Other plots
errorbar(T_dat(:,1),T_dat(:,2),'-o','LineWidth',2,'MarkerFaceColor','auto');
xlim([0 length(T_dat(:,1))+1]);
hold on;
errorbar(C_dat(:,1),C_dat(:,2),'-o','LineWidth',2,'MarkerFaceColor','auto');
legend('TelC','Control');
hold off;ylim([-15 120]);