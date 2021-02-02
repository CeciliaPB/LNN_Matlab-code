
video_list = dir(['*PrL2','*.mp4']);
video_list = {video_list.name}';

% Set the name of the table
A = video_list{1};
totaltab = [genvarname(A(5:end-26)),'_lick'];
tab = table;

%% Licking analysis
% Variables calculated and saved in a table: 
% - tone_Speed:     Speed during the tone;
% - pretone_Speed:  Speed before the tone starts; 
% - light_Speed:    Speed during the light; 
% - prelight_Speed: Speed before the light starts; 
% Change kk to the video you want to analyse

for kk = 1:10

    video = video_list{kk};
    Trackfile = [video(1:8),'_track.mat']; % Set file to save the time variable
    TTLfile   = [video(1:8),'_TTLs.mat']; % Set the file to save
    if kk == 12 || kk == 15
        Trackfile = [video(1:7),'_Test_Track.mat']; % Set file to save the time variable
        TTLfile   = [video(1:7),'_Test_TTLs.mat']; % Set the file to save
    end
    load(TTLfile,'TTL_shock','TTL_light','TTL_sound','delay','t_dig');
    load(Trackfile,'licking_inout','licking_inout2');
    SessionData = table; % Table for each session
    tabname = video(1:end-26);
    shock = TTL_shock - delay;
    sound = TTL_sound - delay;
    light = TTL_light - delay;

%    licking_inout2 = zeros(1,length(licking_inout));
%    for ii = 1:length(licking_inout)-2000
%     A = sum(licking_inout(ii:ii+2000));
%     if A > 500
%       licking_inout2(ii) = 1;
%     elseif A <= 500
%       licking_inout2(ii) = 0;
%     end
%    end

    for ii = 1:length(sound)

        SessionData.Track(ii,:) =[tabname '_s' num2str(ii,'%02.0f')];

        % Set the times of the TTL for each trial
        if exist('sound','var')
            sd = sound(ii)-0.2;
        end

        if exist('light','var')
            li = light(ii)-0.2;
        end

        if exist('shock','var')
            sk = sound(ii)-0.2+18;
        end

        idx = find(isnan(licking_inout2));
        t2 = downsample(t_dig,1000);
%         inout2 = downsample(licking_inout2,1000);
        inout2 = licking_inout2;
        if     isempty(idx) == 1
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

        % Sound-related data ----------------------------------------------
        if exist('sd','var')
            % To extract time the animal is licking during the tone
            toneDur   = (t2>=sd & t2<=(sd+20)); toneDur(1,end) = 0; % hold on; plot((t2),toneDur*0.06,'g'); hold off
            toneIO = inout2 & toneDur;
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

            % To extract time the animal is licking 90s pretone
            pretoneDur   = (t2>=(sd-90) & t2<=sd); pretoneDur(1,end) = 0;
            pretoneIO = (inout2) & pretoneDur;
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

            % To extract time the animal is licking all trial
            allDur = (t2>=(sd-90) & t2<=(sd+20)); allDur(1,end) = 0;
            allIO = inout2 & allDur;
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

            if     isempty(tone_t)== 1
                tone_t = [0 0];
            elseif isempty(tone_d) == 1
                tone_d = [0 0];
            elseif isempty(pretone_t)== 1
                pretone_t = [0 0];
            elseif isempty(all_t)== 1
                all_t = [0 0];
            end

            SessionData.lick_tone(ii,:)    = sum(tone_t(:,2)-tone_t(:,1))/(tone_d(:,2)...
                                            -tone_d(:,1));
            SessionData.lick_pretone(ii,:) = sum(pretone_t(:,2)-pretone_t(:,1))...
                                            /(pretone_d(:,2)-pretone_d(:,1));
            SessionData.lick_all(ii,:)      = sum(all_t(:,2)-all_t(:,1))/(all_d(:,2)-all_d(:,1));
        end

        % Light-related data ----------------------------------------------
        if exist('li','var')
            % To extract time the animal is licking during the light
            lightDur = (t2>=li & t2<=(li+10)); lightDur(1,end) = 0; % hold on; plot((t2),lightDur*0.06,'g'); hold off
            lightIO = (inout2) & lightDur;
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

            if isempty(light_t)== 1
                light_t = [0 0];
            end 

            SessionData.lick_light(ii,:)    = sum(light_t(:,2)-light_t(:,1))...
                                            /(light_d(:,2)-light_d(:,1));
        end

    end

    tab = [tab; SessionData];
    eval([totaltab '= [tab];']);

end

%% Add to Control/ TelC
% My inhibited (TelC) animals are #1, 3, 7, 8. Inhibited = [1, 3, 7, 8];
% My control (Control) animals are #2, 4, 5, 6. Control = [2, 4, 5, 6];
% Change the name of the structure and the animals to analyse accordingly.

parentpath = pwd;
childpath = dir('OF*');
childpath = {childpath.name}';

% Inhibited = [1, 3, 7, 8];
Inhibited = [5,6];
for ii = Inhibited
    cd([parentpath,filesep,childpath{ii}])
    Animal = dir(['OF*', num2str(ii), '*_data.mat']);
    Animal = {Animal.name}';
    load(Animal{1,:});
    
    s = whos('*lick*'); 
    Tablename = s.name;
    Data = evalin('base',[Tablename '(:,:);']);
    
    VarName = Data.Properties.VariableNames;
    
    if ii == Inhibited(1)
        Inhib.(['Lick' VarName{2}]) = Data.(VarName{2});     
        Inhib.(['Lick' VarName{3}]) = Data.(VarName{3}); 
        Inhib.(['Lick' VarName{4}]) = Data.(VarName{4});     
        Inhib.(['Lick' VarName{5}]) = Data.(VarName{5}); 
    elseif ii ~= Inhibited(1)
        Inhib.(['Lick' VarName{2}]) = [Inhib.(['Lick' VarName{2}]),Data.(VarName{2})];     
        Inhib.(['Lick' VarName{3}]) = [Inhib.(['Lick' VarName{3}]),Data.(VarName{3})]; 
        Inhib.(['Lick' VarName{4}]) = [Inhib.(['Lick' VarName{4}]),Data.(VarName{4})];     
        Inhib.(['Lick' VarName{5}]) = [Inhib.(['Lick' VarName{5}]),Data.(VarName{5})]; 
    end
    
    clearvars -except Ctrl Inhib parentpath childpath Control Inhibited
    cd(parentpath)
end

%% Ploting the data

Data = PrL1_lick;
Animal_Data = struct;
nn = 1; % Grouping number, to analyse data in groups of nn trials

Animal_Data.lick_tone     = zeros(length(Data.Track)/nn,1);
Animal_Data.lick_pretone  = zeros(length(Data.Track)/nn,1);
Animal_Data.lick_light    = zeros(length(Data.Track)/nn,1);
Animal_Data.lick_all      = zeros(length(Data.Track)/nn,1);
for ii = 1:nn:length(Data.Track)
    lick_tone     = nanmean(Data.lick_tone(ii:ii+(nn-1),:));
    lick_pretone  = nanmean(Data.lick_pretone(ii:ii+(nn-1),:));
    lick_light    = nanmean(Data.lick_light(ii:ii+(nn-1),:));
    lick_all      = nanmean(Data.lick_all(ii:ii+(nn-1),:));
    
    n = ceil(ii/nn);
    Animal_Data.lick_tone(n,1)     = lick_tone;
    Animal_Data.lick_pretone(n,1)  = lick_pretone;
    Animal_Data.lick_light(n,1)    = lick_light;
    Animal_Data.lick_all(n,1)      = lick_all;
    
end

figname = Data.Track(1,5:end-3);

plot(Animal_Data.lick_tone*100,'b-o','LineWidth',2,'MarkerFaceColor','b');
hold on;
plot(Animal_Data.lick_pretone*100,'b--o','LineWidth',2,'MarkerFaceColor','b');
plot(Animal_Data.lick_light*100,'c-o','LineWidth',2,'MarkerFaceColor','c');
plot(Animal_Data.lick_all*100,'r-o','LineWidth',2,'MarkerFaceColor','r');
hold off;
ylim([-5 105]);
xlim([0 length(Animal_Data.lick_tone)+1]);
legend('Time INplat - tone','Time INplat - pretone','Time INplat - light',...
    'Time INplat - total','Location', 'northeastoutside');
title('Time spent licking');
xlabel('Trial #'); 
ylabel('% Of Time'); 
hold off;
% set(gca,'fontname','arial');  % Set font to arial
% saveas(gca,[figname, '_TimeIN'],'png');

save([figname, 'Data_lick'],'Animal_Data','-append');
