%% VIDEO ANALYSIS: Steps to perform.
% Requires the Computer Vision System & Image Processing Toolboxes
% Remember to comment ALL steps you don't need
sd% The parameters to adjust are (in order of relevance): 1) The DETECTOR;...
% 2) The MASKs; 3) The BLOB.

% Read the video and plot the audio to select the times to analyse
video = 'WIN_20200218_10_14_53_Pro.mp4';
[~,name,~] = fileparts(video);
audio = audioread(video);
info = audioinfo(video);
t = linspace(0,info.Duration,length(audio));
plot (t,audio);

% Save the audio file
% audio = typecast(audio(:,1), 'uint64'); % To later plot this convert to double
% t = typecast(t, 'uint64');
save([name '_audio'],'audio','t');

name = video;
% audio = typecast(audio, 'double'); t = typecast(t, 'double');

%% Filter the audio in case of noisy background sound
fs = 44100;
wn=[10000 12000]/(fs/2);   %bandpass 7495 7505 for sound duration
[b,a]=butter(3,wn);  
f=filter(b,a,audio(:,1));
N = 4096;
freq = linspace(0,fs,N);
F = fft(f,N);
maxFreq = N/8; %~2756 Hz.
% figure; subplot(2,1,1)
figure; plot(t, f)
% subplot(2,1,2)
% plot(freq(1:maxFreq),abs(F(1:maxFreq)));

%% Detect sound start times and set end
% load('A9_T10_audio.mat');
% audio = typecast(audio, 'double'); t = typecast(t, 'double');
info = audioinfo(name);
audio = audioread(name);
t = linspace(0,info.Duration,length(audio));

fs = 44100;
wn=[2500 3500]/(fs/2);   %bandpass 7495 7505 for sound duration
[b,a]=butter(3,wn);  
f=filter(b,a,audio);
% plot(t, f)

fm = mean(f);
fsd = std(f);
ft = t>50;
fth = f>fm+fsd*2;

[peaks,locs] = findpeaks(double(fth(ft)),t(ft),'MinPeakDistance',35);
sound = transpose(locs + 30);

% clearvars -except sound
% hold on;
% plot(sound,peaks*0.005,'.')
% hold off;

%% Select the period to analyse and save the video you will use later 
%  - NOT NEEDED NOW!
period = [853 1033]; 

% time = ([period(1)-2600 period(2)-2600]);
% cut_videolapse(video,[time(1) time(2)]);

cut_videolapse(video,[period(1) period(2)]);
name = [video(1:end-4),'_', num2str(period(1)), '_', num2str(period(2)),'.avi'];

%% TRACKING! 

% Setting the detector and blob properties ---------------------------------
% IMPORTANT DETECTOR: Change: 'NumTrainingFrames' to increase or decrease
% the background period; 'LearningRate' value between 0.01 (if the animal
% moves fast) to 0.000001 (if the animal is lazy or the detection too
% jumpy); 'InitialVariance', increase if there is noise or decrease if no
% detection; 'MinimumBackgroundRatio', increase if there is background
% noise or decrease if there is no detection.
% IMPORTANT BLOB: Modify the BlobArea, min and max. Objects over (max)... 
% or under (min) values will not be considered
detector = vision.ForegroundDetector('NumTrainingFrames',15,... 
    'AdaptLearningRate',1,'LearningRate',0.000001,... % Recomended LR 0.005 to 0.0000001
    'InitialVariance',(80/250)^2,'MinimumBackgroundRatio',0.4);  

blob = vision.BlobAnalysis('CentroidOutputPort',true,'AreaOutputPort',true,...
    'BoundingBoxOutputPort', true,'MinorAxisLengthOutputPort',true,...    
    'MinimumBlobAreaSource', 'Property', 'MinimumBlobArea', 150,...
    'MaximumBlobArea',30000);
   
shapeInserter = vision.ShapeInserter('Shape','Circles','BorderColor',...
    'Custom','CustomBorderColor',[255 255 0],'LineWidth',2); % 'Circles' or 'Rectangles'

% Tracking ----------------------------------------------------------------
% IMPORTANT MASKs: A series of modifications of the detected object...
% Named: mask to mask5 and frMask. See section: For detection of moving
% object modify the [Row Column] values. Some recomendations in the actual
% part.

% video = [name '.mp4']; % name of the video to analyse
video = name;

videoSource = vision.VideoFileReader(video); % Read the video

track = []; % Pre-allocate for speed
inout = []; % Pre-allocate for speed

frameSize = size(step(videoSource));        
videoPlayer = vision.VideoPlayer('Position',[100 100 [frameSize(2),...
    frameSize(1)]+30]); % Setting the dimension of the frame displayed

% Select a ROI, manual entry. 
% Also possible: load a inArea bolean as ROI, I need to write this part.
% First ROI for cage detection
Area = exist('inArea','var'); % If a 'inArea' already exists you don't need... 
% to enter the ROI again, CAREFUL! when changing videos be sure to delete this variable
if Area == 1
elseif Area ~= 1
    [X,Y] = meshgrid(1:frameSize(2),1:frameSize(1));
    frame = step(videoSource);
    figure('name','Area'); imshow(step(videoSource));
    [x,y] = ginput;
    inArea = inpolygon(X(:),Y(:),x,y);
    inArea = double(reshape(inArea,[frameSize(1) frameSize(2)]));
    frame(~inArea) = 0;
    hold on; imshow(frame); hold off;
    close Area;
end

% Second ROI for platform detection. Comment if not used, then uncomment
% OUT variable (see later)
Plat = exist('inPlat','var'); % Idem as for 'inArea'
if Plat == 1    
elseif Plat ~= 1
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

% % Uncomment only to save video
% videoFWriter = VideoWriter('Tracking.mp4','MPEG-4'); % Save mp4
% videoFWriter.FrameRate =  videoSource.info.VideoFrameRate;
% videoFWriter.Quality = 100;
% open(videoFWriter);
% vw = figure('name','vw','units','points','position',[20,20,600,585]);
% subplot(4,1,[1 2 3]);
% subplot(4,1,4); set(vw,'DefaultLegendAutoUpdate','off')
%         sk = shock(1)-0.2; % Shock time up to 1st decimal
%         sd = sound(1)-0.2;
%         li = light(1)-0.2;
%         area(([sk sk+2]),([1.3 1.3]),'FaceAlpha',0.5,'FaceColor',...
%         [0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840]);
%         hold on;
%         area(([li li+10]),([1.2 1.2]),'FaceAlpha',0.5,'FaceColor',...
%         [0.9290 0.6940 0.1250],'EdgeColor',[0.9290 0.6940 0.1250]);
%         area(([sd sd+20]),([1.1 1.1]),'FaceAlpha',0.5,'FaceColor',...
%         [0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410]);
%         hold on;
%         sk = shock(2)-0.2; % Shock time up to 1st decimal
%         sd = sound(2)-0.2;
%         li = light(2)-0.2;
%         area(([sk sk+2]),([1.3 1.3]),'FaceAlpha',0.5,'FaceColor',...
%         [0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840]);
%         area(([li li+10]),([1.2 1.2]),'FaceAlpha',0.5,'FaceColor',...
%         [0.9290 0.6940 0.1250],'EdgeColor',[0.9290 0.6940 0.1250]);
%         area(([sd sd+20]),([1.1 1.1]),'FaceAlpha',0.5,'FaceColor',...
%         [0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410]);
%         xlim([0 sd+22]);
%         hold off;
%         legend('Shock','Light','Sound','Location','northwest');
%         yticks([0 1])
%         yticklabels({'OUT','IN'})
%         ylim([0 1.5]);
%         ii = 0;

while ~isDone(videoSource)
     fr  = videoSource();     % Read frame
%      fr = imadjust(fr,[],[.2 .9 .9; 1 1 1],2); % Ccolor detection. Enhances blacks
     
%      % For detection of LEDs
%      % CHANGE BLOB:'MinimumBlobArea',~50; DETECTOR:'MinimumBackgroundRatio',0.9
%      mask = logical(sum((fr>0.9),3)); % Binary detection, get px over a thr (0.9)
%      frMask = repmat(mask & inArea,[1 1 1]);
%      [xy,r]   = blob(frMask);
%      track = [track; xy]; % Saves the tracking points, to plot them later.
% %      circ = [xy,r];
%      out    = insertMarker(fr,track,'o','color',[255 255 0],'size',1);
% %      out    = shapeInserter(fr,circ); % Circle on detected object
% %      bbox   = blob(frMask); % For rectangle shapes: vision.BlobAnalysis,...
% %           select 'BoundingBoxOutputPort', true
% %      out    = shapeInserter(fr,bbox); % Rectagle on detected object
%      step(videoPlayer, out);
%       
     % For detection of moving object. Tips for modification added
     fr2 = imadjust(fr,[],[],1); % Enhance gamma to dark (KINGA)
     mask = detector(fr2); % Binary detection, using vision.ForegroundDetector
     mask1 = repmat(mask & inArea,[1 1 1]);
     mask2 = imopen(mask1, strel('rectangle', [10 10]));
     % MASK2 If Noisy detection: [20 10]; If No detection [3 3]. Very noisy
     % detections: change the DETECTOR: MinimumBackgroundRatio up to 0.9 or
     % increase 'InitialVariance' to (90/250)^2.
     mask3 = imclose(imdilate(edge(mask2,'canny'), strel('rectangle',[5 5])),...
         strel('rectangle', [10 15])); 
     % MASK3 imdilate to make objects bigger [10 10]; imclose to close
     % neighbor pixels: [15 20] for bigger closure.
     mask4 = imfill(mask3, 'holes');
     mask5 = bwareaopen(mask4, 50);
     % MASK5 Remove small objects, adjust to your detection.
     frMask = imerode(mask5, strel('rectangle', [5 3]));
     % frMASK Erosion of the objects, bigger erosion [20 30], smaller erosion [5 3].
     [a,xy,bbox,r]   = blob(frMask);
    
     A = a == max(a); % Selects the biggest blob (usually the animal)
     if isempty(xy)==1 
         xy2 = xy;
         r2 = r;
         nan = [NaN, NaN];
         track = [track; nan]; % Saves the tracking points.
         circ = [xy2,r2];
     elseif isempty(xy)==0  
         xy2 = xy(A,:);
         r2 = r(A,:);
         track = [track; xy2]; % Saves the tracking points, to plot them later.
         circ = [xy2,r2];
     end
     
     
     % Comment this section if there is not second ROI
     if isempty(xy)==0
         % Logical to place the coordinates of the detected object
         position = zeros(frameSize(1),frameSize(2));
         position(round(xy2(1,2)),round(xy2(1,1)))= 1;
         INplat = repmat(position & inPlat,[1 1 1]);
         position(INplat == 0) = [];

         if isempty(position)==0
             out = insertMarker(fr,xy2,'+','color',[255 0 0],'size',5);
             inout = [inout; 1]; % The mice is IN the platform
         elseif isempty(position)==1
             out = shapeInserter(fr,circ); % Circle on detected object
             inout = [inout; 0]; % The mice is OUT the platform
         end

     elseif isempty(xy)==1
         out    = shapeInserter(fr,circ);
         inout = [inout; NaN];
     end
      
% %      % Only to save film
%      vw = figure('name','vw','units','points','position',[20,20,600,585]);
% subplot(4,1,[1 2 3]);
% subplot(4,1,4); set(vw,'DefaultLegendAutoUpdate','off')
%         sk = shock(1)-0.2; % Shock time up to 1st decimal
%         sd = sound(1)-0.2;
%         li = light(1)-0.2;
%         area(([sk sk+2]),([1.3 1.3]),'FaceAlpha',0.5,'FaceColor',...
%         [0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840]);
%         hold on;
%         area(([li li+10]),([1.2 1.2]),'FaceAlpha',0.5,'FaceColor',...
%         [0.9290 0.6940 0.1250],'EdgeColor',[0.9290 0.6940 0.1250]);
%         area(([sd sd+20]),([1.1 1.1]),'FaceAlpha',0.5,'FaceColor',...
%         [0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410]);
%         hold on;
%         sk = shock(2)-0.2; % Shock time up to 1st decimal
%         sd = sound(2)-0.2;
%         li = light(2)-0.2;
%         area(([sk sk+2]),([1.3 1.3]),'FaceAlpha',0.5,'FaceColor',...
%         [0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840]);
%         area(([li li+10]),([1.2 1.2]),'FaceAlpha',0.5,'FaceColor',...
%         [0.9290 0.6940 0.1250],'EdgeColor',[0.9290 0.6940 0.1250]);
%         area(([sd sd+20]),([1.1 1.1]),'FaceAlpha',0.5,'FaceColor',...
%         [0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410]);
%         xlim([0 sd+22]);
%         hold off;
%         legend('Shock','Light','Sound','Location','northwest');
%         yticks([0 1])
%         yticklabels({'OUT','IN'})
%         ylim([0 1.5]);
% 
%      figure(vw); hold on; 
%      subplot(4,1,[1 2 3]); imshow(fr);
%      imshow(out);
%      hold off;
%      hold on;
%      plot(track(:,1),track(:,2),'g','LineWidth',1);
%      subplot(4,1,4);  
%      t3 = linspace(period(1),...
%          period(1)+(length(inout)/videoSource.info.VideoFrameRate),...
%          length(inout));
%      plot(t3,inout,'k');
%      hold off;
%      ii = ii + 1;
%      saveas(vw, genvarname(['MOV',num2str(ii)]), 'png');
%      close(vw) 
%      videoFrame = getframe(vw);
% %      set(vw,'PaperType','a3');
% %      videoFrame = print ('-RGBImage','-r300');
%      videoFWriter(videoFrame.cdata);

% % Uncomment the one you need in case there is no second ROI     
%       out = insertMarker(fr,xy2,'+','color',[255 255 0],'size',5);
% %      out    = shapeInserter(fr,circ); % Circle on detected object
% %      out    = shapeInserter(fr,bbox); % Rectagle on detected object

     step(videoPlayer, out);

end
% close(videoFWriter); % To write the video, uncommet to save the video
release(videoPlayer);
clearvars -except name audio t detector inout t2 track trigger video_list inArea

%% Plot the audio and position INOUT platform of the animal

 audio = typecast(audio(:,1), 'double'); t = typecast(t, 'double');
% t2 = linspace(period(1)-(5/videoSource.info.VideoFrameRate),period(2),length(inout));
% t2 = linspace(period(1),period(2),length(inout));
figure; plot (t,audio,'b');
% xlim([period(1) period(2)]);
hold on;
% plot((t2),inout*0.05,'r');
t2 = linspace(min(t),max(t),length(inout));
plot((t2),inout*0.05,'r'); 
ylim([-0.05 0.06]);
hold off;

%% Saving
FIG = input('Name?:','s');

% save([FIG '.mat'],'detector','inout','period','t2','track');
save([FIG '.mat'],'detector','inout','t2','track','trigger'); 
%   saveas(gcf, FIG, 'fig');
saveas(gcf, FIG, 'png');

%% Option B: When analysing whole-Test video 
for ii = 1:length(trigger)
fin = trigger(ii);
 
period = [ceil(fin)-178 ceil(fin)+2];
% figure; plot (t,audio2,'b');
xlim([period(1) period(2)]);
% hold on;
% plot((t2),inout*0.05,'r'); 
% hold off;

FIG2 = [FIG,'_',num2str(period(1)),'_',num2str(period(2))];
saveas(gcf, FIG2, 'png');
end

%% Matlab video player, only plays the frames
% videoSource = vision.VideoFileReader(name);
% videoPlayer = vision.VideoPlayer(); 
% while ~isDone(videoSource)
% frame  = videoSource();
% videoPlayer(frame);
% end

%% Time spent in platform
% Variable 'trigger' is a list of the end time of the sound. Change table
% name according to your Animal/ Test.

A9_T10 = table;
tab = whos('A9_T*');
A9_T10.Track =[[tab.name '_s11'];[tab.name '_s12'];[tab.name '_s13'];...
    [tab.name '_s21'];[tab.name '_s22'];[tab.name '_s23'];...
    [tab.name '_s31'];[tab.name '_s32'];[tab.name '_s33']];

for ii = 1:length(trigger)
    fin = trigger(ii); % Last tone time

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
    toneDur = (t2>=(fin-30) & t2<=fin); toneDur(1,end) = 0; % hold on; plot((t2),toneDur*0.06,'g'); hold off
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
    
    % And the latency to enter the platform
    if isempty(B)
        min_lat = 30;
        max_lat = 30;
    elseif ~isempty(B)
        min_lat = tone_t(1,1) - tone_d(1,1);
        max_lat = tone_t(end,1) - tone_d(1,1);
    end
    
    % To extract time the animal is inout the platform 2min pretone
    timeDur = (t2>=(fin-150) & t2<=fin-30); timeDur(1,end) = 0;
    preIO = transpose(inout2) & timeDur;
    preIO(1,1) = 0;

    A = ischange(double(preIO),'MaxNumChanges',length(preIO));
    B = t2(A); 
    if mod(length(B),2)==0 
    elseif mod(length(B),2)~=0
        B = [B,t2(end)];
    end
    pretone_t = transpose(B(:,1:2:end));  % odd matrix
    pretone_t(:,2) = transpose(B(:,2:2:end)); % even matrix

    A = ischange(double(timeDur),'MaxNumChanges',length(timeDur));
    pretone_d = t2(A);

    % To extract time the animal is inout the platform all
    allDur = (t2>=(fin-150) & t2<=fin); allDur(1,end) = 0;
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
    shockDur = (t2>=(fin-2) & t2<=fin); % hold on; plot((t2),toneDur*0.06,'g'); hold off
    shockIO = transpose(inout2) & shockDur;
    S = find(shockDur == 1);
    if (sum(abs(shockIO(S)-1))) >= 2
         shock = 1;
    else 
        shock = 0;
    end 

    % Calculations
    if isempty(tone_t)== 1
        tone_t = [0 0];
    elseif isempty(pretone_t)== 1
        pretone_t = [0 0];
    elseif isempty(all_t)== 1
        all_t = [0 0];
    end 

    A9_T10.tIN_tone(ii,:) = 100*sum(tone_t(:,2)-tone_t(:,1))/(tone_d(:,2)-tone_d(:,1));
    A9_T10.tIN_pretone(ii,:) = 100*sum(pretone_t(:,2)-pretone_t(:,1))/(pretone_d(:,2)-pretone_d(:,1));
    A9_T10.tIN_all(ii,:) = 100*sum(all_t(:,2)-all_t(:,1))/(all_d(:,2)-all_d(:,1));
    A9_T10.tIN_shock(ii,:) = shock;
    A9_T10.min_lat(ii,:) = min_lat;
    A9_T10.max_lat(ii,:) = max_lat;
        
end

 save([tab.name '_analysis.mat'],tab.name);

%% Now the plot
% Session mean

newtab = load('BA2b_table.mat');
newtab = struct2cell(newtab);
Data = newtab{1,1};
Animal_Data = struct;
nn = 7;

Animal_Data.tone = zeros(length(Data.Track)/nn,1);
Animal_Data.pretone = zeros(length(Data.Track)/nn,1);
Animal_Data.all = zeros(length(Data.Track)/nn,1);
Animal_Data.shock= zeros(length(Data.Track)/nn,1);
for ii = 1:nn:length(Data.Track)
    tone = nanmean(Data.tIN_tone(ii:ii+(nn-1),:));
    pretone = nanmean(Data.tIN_pretone(ii:ii+(nn-1),:));
    all = nanmean(Data.tIN_all(ii:ii+(nn-1),:));
    shock = max(Data.tIN_shock(ii:ii+(nn-1),:));
    n = ceil(ii/nn);
    Animal_Data.tone(n,1) = tone;
    Animal_Data.pretone(n,1) = pretone;
    Animal_Data.all(n,1) = all;
    Animal_Data.shock(n,1) = shock;
    
    min_lat = nanmean(Data.min_lat(ii:ii+(nn-1),:));
    max_lat = nanmean(Data.max_lat(ii:ii+(nn-1),:));
    n = ceil(ii/nn);
    Animal_Data.min(n,1) = min_lat;
    Animal_Data.max(n,1) = max_lat;
end

plot(Animal_Data.tone,'-o','LineWidth',2,'MarkerFaceColor','auto');
hold on;
plot(Animal_Data.pretone,'-o','LineWidth',2,'MarkerFaceColor','auto');
plot(Animal_Data.all,'-o','LineWidth',2,'MarkerFaceColor','auto');
hold off;
ylim([-5 105]);
xlim([0 length(Animal_Data.tone)+1]);

shock = Animal_Data.shock;
shock (shock == 0)= NaN;
hold on;
plot(shock,'o','MarkerEdgeColor',[1 0.41 0.16],'Markersize',10,'MarkerFaceColor',[1 0.41 0.16]);
legend('Time in platform','Time IN pretone','Time IN all','Shock');
hold off;

figure; plot(Animal_Data.min,'-o','LineWidth',2,'MarkerFaceColor','auto');
hold on; plot(Animal_Data.max,'-o','LineWidth',2,'MarkerFaceColor','auto'); 
xlim([0 length(Animal_Data.max)+1]);
legend('Min latency','Max latency');
hold off;

save('A2_lat.mat','Animal_Data');

%% Plot by Group
Data = Control;

rControl.tone = zeros(15,2);
rControl.shock= zeros(15,1);
rControl.shock100 = zeros(15,1);
rControl.min = zeros(15,2);
rControl.max = zeros(15,2);

for ii = 1:3:90
    tone = nanmean(nanmean(Data.tone(ii:ii+2,:)),2);
    shock = max(max(Data.shock(ii:ii+2,:)));
    shock100 = nanmean(sum(Data.shock(ii:ii+2,:)/3*100),2);
    min_lat = nanmean(nanmean(Data.min_lat(ii:ii+2,:)),2);
    max_lat = nanmean(nanmean(Data.max_lat(ii:ii+2,:)),2);
    
    sd_tone = nanmean(nanstd(Data.tone(ii:ii+2,:)),2);
    sd_min_lat = nanmean(nanstd(Data.min_lat(ii:ii+2,:)),2);
    sd_max_lat = nanmean(nanstd(Data.max_lat(ii:ii+2,:)),2);
    
    n = ceil(ii/3);
    rControl.tone(n,1) = tone;
    rControl.tone(n,2) = sd_tone;
    rControl.shock(n,1) = shock;
    rControl.shock100(n,1) = shock100;
    rControl.min(n,1) = min_lat;
    rControl.min(n,2) = sd_min_lat;
    rControl.max(n,1) = max_lat;
    rControl.max(n,2) = sd_max_lat;
end

plot(rControl.tone(:,1),'-o','LineWidth',2,'MarkerFaceColor','auto');
ylim([-5 105]);
xlim([0 length(rControl.tone(:,1))+1]);

shock = rControl.shock;
shock (shock == 0)= NaN;
hold on;
plot(shock,'o','Markersize',10,'MarkerFaceColor',[1 0.41 0.16]);
legend('Time in platform','Shock');
hold off;

figure; plot(rControl.min(:,1),'-o','LineWidth',2,'MarkerFaceColor','auto');
hold on; 
plot(rControl.max(:,1),'-o','LineWidth',2,'MarkerFaceColor','auto');
legend('Min latency','Max latency');
hold off;

%% T-test between TelC and Controls by Day
DataT = TelC.licks_dur;
DataC = Control.licks_dur;
name = 'LicksDur_dat';
k = 9; % Set to 3 for session, to 9 for day analysis

p_val = zeros(5,1);
C_dat = zeros(5,2);
T_dat = zeros(5,2);

for ii = 1:k:90
    
mT = nanmean(DataT(ii:ii+(k-1),:),2);
mC = nanmean(DataC(ii:ii+(k-1),:),2);
[h,p] = ttest2(nanmean(DataT(ii:ii+(k-1),:),2),nanmean(DataC(ii:ii+(k-1),:),2));

n = ceil(ii/k);
p_val(n,1) = p;
T_dat(n,1) = nanmean(mT);
T_dat(n,2) = nanstd(mT);
C_dat(n,1) = nanmean(mC);
C_dat(n,2) = nanstd(mC);

end

pval = double(p_val<=0.05)*120;
pval(pval == 0)= NaN;
errorbar(T_dat(:,1),T_dat(:,2),'-o','LineWidth',2,'MarkerFaceColor','auto');
xlim([0 length(T_dat(:,1))+1]);
hold on;
errorbar(C_dat(:,1),C_dat(:,2),'-o','LineWidth',2,'MarkerFaceColor','auto');
plot(pval(:,1),'k*');
legend('TelC','Control','p-val < 0.05');
ylim([-15 130]);
hold off;

save(name,'T_dat','C_dat','p_val');
