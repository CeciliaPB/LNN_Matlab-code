function wf = waveform_analysis_spike2(varargin)

% This function calculates the features of the waveforms and/ or plots the
% waveforms of the channel where the spike is detected.
% 
% INPUTS: 
%   - A window will open to choose the .mat file with the mean waveform of
%   the neuron. ex. 'GR1_1.mat' 
%   Varargin 
%   - 'plot': When you want to plot the waveform. The default plot is the
%   mean waveform.  
%   - 'features': generates the structure wf, that contains all the
%   calculated features. 
%   - 's': saves the figure in the specified file format (jpg, png, svg,
%   etc.). See saveas help for more info. Default NOT save.
%
% OUTPUTS: 
%   - wf: a structure with the calculated features. For more information
%   about these features click to see <a href="matlab:A =
%   imread('L:\Cecilia\WFparams.tif'); imshow(A);">Features image</a>. 
%   - without output arguments, the function plots the mean waveform
% 
% Examples: 
% wf = waveform_analysis('features');
% wf = waveform_analysis('features','plot','s','jpg');
% waveform_analysis(group,neuron); 
% 
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2019
% Laboratory of Network Neurophysiology
% Institute of Experimantal Medicine, Hungary.
%
% MATLAB toolboxes: Signal Processing Toolbox.
% -------------------------------------------------------------------------

clearvars -except varargin 
% 1_ Defining params ------------------------------------------------------
s = 0;

% Params, read from varargin
if nargin
  for ii=1:size(varargin,2)
    switch varargin{ii}
      case 's'
        s = 1;
        type = varargin{ii+1};
      otherwise  
    end
  end
else
end

% 2_ load the .mat file ---------------------------------------------------
[name, ~] = uigetfile('*.mat*', 'Select a Mat File', 'MultiSelect', 'off');
    if (name == 0)
      return;
    end
WF = load(name);
WFfileds = fieldnames(WF);
WF = WF.(WFfileds{1,1});

WF.wflength = WF.length;
WF = rmfield(WF,{'length','title'});

WFfileds = fieldnames(WF);
for kk=1:length(WFfileds)
    eval([WFfileds{kk} '=WF.' WFfileds{kk}]);
end

% 3_ generate the vars used to calculate ----------------------------------
mUP = values(:,:); %#ok<IDISVAR,NODEF> % mean wf of the neuron
mDOWN = -mUP; % Inverted values
time = linspace(start, start + (interval * wflength),wflength);

Type = imread('L:\Cecilia\WFparams_Types.tif');
figure('name','WaveformType','Position',[100 200 [838, 600]]);
subplot(2,1,1);
plot(time,values,'b','LineWidth',2,'Color',[0.8500 0.3250 0.0980]);
title('Your waveform');
subplot(2,1,2,'Position',[.0 -.25 [1, 1]]);
imshow(Type);
disp('Please select the waveform type,');
TypeID = input('[1/ 2]:');
close WaveformType

if logical(max(strcmp(varargin, 'features')) == 1) && ...
    logical(max(strcmp(varargin, 'plot')) == 1)
    
    % Calculate features using the mean waveform, depends on waveform type
    pk = zeros(1,2);
    loc = zeros(1,2);
    pr = zeros(2,2);
    switch TypeID
        case 1
            [pk2(1,1),pk2(2,1),w2,~] = findpeaks(mUP,time,'SortStr',...
                'descend','NPeaks',1);
            [p,l,~,r] = findpeaks(mDOWN,time,'MinPeakHeight',0.01);
        case 2
            [pk2(1,1),pk2(2,1),w2,~] = findpeaks(mDOWN,time,'SortStr',...
                'descend','NPeaks',1);
            [p,l,~,r] = findpeaks(mUP,time,'MinPeakHeight',0.01);
            pk2(1,1) = -pk2(1,1);
            p = -p;
    end
    
    % We want these variables being a 1 row vector.
    if size(p,1) > size(p,2)
        p = transpose(p);
    end
    if size(l,1) > size(l,2)
        l = transpose(l);
    end
    if size(r,1) > size(r,2)
        r = transpose(r);
    end
    
    if size(p,2) == 1 || size(l,2) == 1 || size(r,2) == 1
        pk(:,1) = p;
        loc(:,1) = l;
        pr(:,1)= r;
    elseif size(p,2) >= 2
        A = find(l < pk2(2,1));
        B = find(l > pk2(2,1));
        try
            pk(1,1) = p(A(end));
        catch
            pk(1,1) = 0;
        end
        pk(1,2)  = p(B(1));
        try
            loc(1,1) = l(A(end));
        catch
            loc(1,1) = 0;
        end
        loc(1,2) = l(B(1));
        try
            pr(1,1) = r(A(end));
        catch
            pr(1,1) = 0;
        end        
        pr(1,2)  = r(B(1));
    end
    
    % Calculate values
    F = mDOWN>0;
    I(1,1) = find(F,1,'first');
    I(1,2) = find(F,1,'last');
    wf.a = pr(1,1);
    wf.b = pr(1,2);
    wf.pk1(1,1) = -pk(1,1);
    wf.pk1(1,2) = loc(1,1);
    wf.pk2 = transpose(pk2);
    wf.pk3(1,1) = -pk(1,2);
    wf.pk3(1,2) = loc(1,2);
    wf.c = wf.pk2(1,2) - wf.pk1(1,2);
    wf.d = wf.pk3(1,2) - wf.pk2(1,2);
    wf.e = abs(wf.pk2(1,1)) + abs(wf.pk1(1,1));
    wf.f = abs(wf.pk3(1,1)) + abs(wf.pk2(1,1));
    wf.g = w2;
    wf.first = wf.pk2(1,2) - I(1,1);
    wf.last  = I(1,2) - wf.pk2(1,2);
    wf.mWF   = mUP;
    wf.time  = time;
    
    % plot
    plot(time*1000,mUP);
    xlabel('Time (ms)');
    ylabel('Amplitude (mV)');
    title('Waveform');
    set(gca,'fontname','arial');
    hold on;
    plot(wf.pk1(1,2)*1000,wf.pk1(1,1),'o','Color',[0.3010 0.7450 0.9330])
    plot(wf.pk2(1,2)*1000,wf.pk2(1,1),'o','Color',[0.6350 0.0780 0.1840])
    plot(wf.pk3(1,2)*1000,wf.pk3(1,1),'o','Color',[0.8500 0.3250 0.0980])
    legend('Waveform','pk1','pk2','pk3','Location','northeast');

elseif max(strcmp(varargin, 'plot')) == 1
    plot(time*1000,mUP);
    xlabel('Time (ms)');
    ylabel('Amplitude (mV)');
    title('Waveform');
    set(gca,'fontname','arial');
    
elseif max(strcmp(varargin, 'features')) == 1
    
    % Calculate features using the mean waveform, depends on waveform type
    pk = zeros(1,2);
    loc = zeros(1,2);
    pr = zeros(2,2);
    switch TypeID
        case 1
            [pk2(1,1),pk2(2,1),w2,~] = findpeaks(mUP,time,'SortStr',...
                'descend','NPeaks',1);
            [p,l,~,r] = findpeaks(mDOWN,time,'MinPeakHeight',0.01);
        case 2
            [pk2(1,1),pk2(2,1),w2,~] = findpeaks(mDOWN,time,'SortStr',...
                'descend','NPeaks',1);
            [p,l,~,r] = findpeaks(mUP,time,'MinPeakHeight',0.01);
            pk2(1,1) = -pk2(1,1);
            p = -p;
    end
    
    % We want these variables being a 1 row vector.
    if size(p,1) > size(p,2)
        p = transpose(p);
    end
    if size(l,1) > size(l,2)
        l = transpose(l);
    end
    if size(r,1) > size(r,2)
        r = transpose(r);
    end
    
    if size(p,2) == 1 || size(l,2) == 1 || size(r,2) == 1
        pk(:,1) = p;
        loc(:,1) = l;
        pr(:,1)= r;
    elseif size(p,2) >= 2
        A = find(l < pk2(2,1));
        B = find(l > pk2(2,1));
        try
            pk(1,1) = p(A(end));
        catch
            pk(1,1) = 0;
        end
        pk(1,2)  = p(B(1));
        try
            loc(1,1) = l(A(end));
        catch
            loc(1,1) = 0;
        end
        loc(1,2) = l(B(1));
        try
            pr(1,1) = r(A(end));
        catch
            pr(1,1) = 0;
        end
        pr(1,2)  = r(B(1));
    end
    
    % Calculate values
    F = mDOWN>0;
    I(1,1) = find(F,1,'first');
    I(1,2) = find(F,1,'last');
    wf.a = pr(1,1);
    wf.b = pr(1,2);
    wf.pk1(1,1) = -pk(1,1);
    wf.pk1(1,2) = loc(1,1);
    wf.pk2 = transpose(pk2);
    wf.pk3(1,1) = -pk(1,2);
    wf.pk3(1,2) = loc(1,2);
    wf.c = wf.pk2(1,2) - wf.pk1(1,2);
    wf.d = wf.pk3(1,2) - wf.pk2(1,2);
    wf.e = abs(wf.pk2(1,1)) + abs(wf.pk1(1,1));
    wf.f = abs(wf.pk3(1,1)) + abs(wf.pk2(1,1));
    wf.g = w2;
    wf.first = wf.pk2(1,2) - I(1,1);
    wf.last  = I(1,2) - wf.pk2(1,2);
    wf.mWF   = mUP;
    wf.time  = time;
        
else % If there are no varargins plot mean waveform
    plot(time*1000,mUP);
    xlabel('Time (ms)');
    ylabel('Amplitude (mV)');
    title('Waveform');
    set(gca,'fontname','arial');
end

disp(['Save as ' name '?']);
saveID = input('[Y/ N]:','s');
switch saveID
    case {'Y', 'y',''}
        save(name(1:end-4), 'wf', '-append');
        if s == 1
            saveas(gcf, genvarname([name(1:end-4),'_wf']), type);
        elseif s ~= 1
        end 
    case {'N', 'n'}
        name2 = input('Choose file name:','s');
        save(name2, 'wf', '-append');
        if s == 1
            saveas(gcf, genvarname([name2(1:end-4),'_wf']), type);
        elseif s ~= 1
        end 
end

end
