function Track_detection

    clear name
    % Look for the video to analyse
    [name, ~] = uigetfile('*.*', 'Select a Video File', 'MultiSelect', 'off');
        if (name == 0)
          return;
        end

    % Choose params
    disp('Use default params?')
    disp('- Nº ROI:             1')
    disp('- Nº Training frames: 1000')
    disp('- Learning rate:      0.000001')
    disp('- Initial Variance:   (80/250)^2')
    default = input('[Y/ N]?:','s');
    
    switch default
        case {'Y', 'y'}
            nROI = 1;
            nTF  = 1000;
            LR   = 0.000001;
            iVar = (80/250)^2;
        case {'N', 'n'}
            nROI = input('Number of ROIs (up to 2)?:');
            nTF = input('Number of Training frames?:');
            LR = input('Learning rate?:');
            iVar = input('Initial Variance?:');
    end
    
    % Setting the detector and blob properties ----------------------------
    % IMPORTANT DETECTOR: Change: 'NumTrainingFrames' to increase or
    % decrease the background period; 'LearningRate' value between 0.01 (if
    % the animal moves fast) to 0.000001 (if the animal is lazy or the
    % detection too jumpy); 'InitialVariance', increase if there is noise
    % or decrease if no detection; 'MinimumBackgroundRatio', increase if
    % there is background noise or decrease if there is no detection.
    % IMPORTANT BLOB: Modify the BlobArea, min and max. Objects over
    % (max)... or under (min) values will not be considered
    detector = vision.ForegroundDetector('NumTrainingFrames',nTF,... 
        'AdaptLearningRate',1,'LearningRate',LR,... % Recomended LR 0.005 to 0.0000001
        'InitialVariance',iVar,'MinimumBackgroundRatio',0.4);  

    blob = vision.BlobAnalysis('CentroidOutputPort',true,'AreaOutputPort',...
        true,'BoundingBoxOutputPort', true,'MinorAxisLengthOutputPort',true,...    
        'MinimumBlobAreaSource', 'Property', 'MinimumBlobArea', 50,...
        'MaximumBlobArea',25000);

    shapeInserter = vision.ShapeInserter('Shape','Circles','BorderColor',...
        'Custom','CustomBorderColor',[255 255 0],'LineWidth',2); % 'Circles' or 'Rectangles'

    % Tracking ------------------------------------------------------------
    % IMPORTANT MASKs: A series of modifications of the detected object...
    % Named: mask to mask5 and frMask. See section: For detection of moving
    % object modify the [Row Column] values. Some recomendations in the
    % actual part.

    video = name;

    videoSource = vision.VideoFileReader(video); % Read the video

    track = []; % Pre-allocate for speed
    inout = []; % Pre-allocate for speed

    frameSize = size(step(videoSource));        
    videoPlayer = vision.VideoPlayer('Position',[100 100 [frameSize(2),...
        frameSize(1)]+30]); % Setting the dimension of the frame displayed

    % Select a ROI, manual entry. 
    
    % First ROI for cage detection
    try
        inArea= evalin('base','inArea');
    catch
    end
    Area = exist('inArea','var'); % If a 'inArea' already exists you don't...  
    % need to enter the ROI again.
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
    
    if nROI == 2
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
    end

    while ~isDone(videoSource)
         fr  = videoSource();     % Read frame
    %      fr2 = imadjust(fr,[],[.2 .9 .9; 1 1 1],2); % Ccolor detection. Enhances blacks
     
         % For detection of moving object. Tips for modification added
         fr = imadjust(fr,[],[],1.3); % Enhance gamma to dark (KINGA)
         mask = detector(fr); % Binary detection, using vision.ForegroundDetector
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
         [a,xy,~,r]   = blob(frMask);

         A = a == max(a); % Selects the biggest blob (usually the animal)
         if isempty(xy)==1 
             xy2 = xy;
             r2 = r;
             nan = [NaN, NaN];
             track = [track; nan]; %#ok<*AGROW> % Saves the tracking points.
             circ = [xy2,r2];
         elseif isempty(xy)==0  
             xy2 = xy(A,:);
             r2 = r(A,:);
             track = [track; xy2]; % Saves the tracking points, to plot them later.
             circ = [xy2,r2];
         end

        if nROI == 2
            % Comment this section if there is not second ROI
            if isempty(xy)==0
                % Logical to place the coordinates of the detected object
                position = zeros(frameSize(1),frameSize(2));
                position(round(xy2(1,2)),round(xy2(1,1)))= 1;
                INplat = repmat(position & inPlat,[1 1 1]);
                position(INplat == 0) = [];

                if isempty(position)==0
                    out    = shapeInserter(fr,circ); % Circle on detected object
                    inout = [inout; 1]; % The mice is IN the platform
                elseif isempty(position)==1
                    out = insertMarker(fr,xy2,'+','color',[255 0 0],'size',5);
                    inout = [inout; 0]; % The mice is OUT the platform
                end

            elseif isempty(xy)==1
                out = insertMarker(fr,xy2,'+','color',[255 0 0],'size',5);
                inout = [inout; NaN];
            end
        end
        
        if nROI == 1
            out = insertMarker(fr,xy2,'+','color',[255 0 0],'size',5);
        end
        
    step(videoPlayer, out);

    end
    
    release(videoPlayer);
    clearvars -except name time detector inout track inArea nROI

    info = audioinfo(name);
    time = linspace(0,info.Duration,length(track));

    % Saving
    Savefile = input('Saving Name?:','s');
    
    if nROI == 1
        save([Savefile '.mat'],'detector','time','track', 'inArea');
    elseif nROI == 2
        save([Savefile '.mat'],'detector','inout','time','track', 'inArea'); 
    end
    
end