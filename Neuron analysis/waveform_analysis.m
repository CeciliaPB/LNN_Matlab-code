function wf = waveform_analysis(group,neuron,varargin)

% This function calculates the features of the waveforms and/ or plots the
% waveforms of the channel where the spike is detected.
% 
% INPUTS: 
%   - group: .mat file of the group where the neuron is. ex. 'GR1.mat'
%   - neuron: .mat file with thw timestamps of the neuron. ex. 'GR1_1.mat'
%   Varargin
%   - 'plot': When you want to plot the waveform. The default plot is the
%   mean waveform. You can choose to plot only the mean waveform ('mean') or
%   all the waveforms with overlaped mean waveform ('all').
%   - 'features': generates the structure wf, that contains all the
%   calculated features.
%   - 's': saves the figure in the specified file format (jpg, png,
%   svg, etc.). See saveas help for more info. Default NOT save.
%
% OUTPUTS: 
%   - wf: a structure with the calculated features. For more information
%   about these features click to see <a href="matlab:A =
%   imread('L:\Cecilia\WFparams.tif'); imshow(A);">Features image</a>. 
%   - without output arguments, the function plots the mean waveform
% 
% Examples: 
% wf = waveform_analysis(group,neuron,'features');
% wf = waveform_analysis(group,neuron,'features','plot','mean','s','jpg');
% waveform_analysis(group,neuron); 
% 
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2019
% Laboratory of Network Neurophysiology
% Institute of Experimantal Medicine, Hungary.
%
% MATLAB toolboxes: Signal Processing Toolbox.
% -------------------------------------------------------------------------

% Default params
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

% Load variables unless they are already in the Workspace
if ~exist('TimeStamps','var') == 1
    load(group,'TimeStamps');
end
if ~exist('WaveForms','var') == 1
    load(group,'WaveForms');
end
if ~exist('TS','var') == 1
    load(neuron,'TS');
end

A = TS/10000;
C = zeros(length(TimeStamps),1);

for kk = 1:length(A)
    B = TimeStamps==A(kk,1);
    C(B) = 1;    
end

for n = 1:size(WaveForms,2)
D = WaveForms(:,n,:);
E = D(logical(C),:,:);
F = mean(E(:,:));
F2(n,:) = F; 
end

H = [min(F2,[],2),max(F2,[],2)];
H2 = H(:,2)-H(:,1);
H3 = max(H2);
G = H2 == H3;
D = WaveForms(:,G,:); % Select the wave from G channel
E = D(logical(C),:,:);
allUP = E(:,:); % All wf of the neuron
allDOWN = -allUP; % Inverted values
mUP = mean(allUP); % Mean wf
mDOWN = -mUP; % Inverted values

if logical(max(strcmp(varargin, 'features')) == 1) && logical(max(strcmp(varargin, 'plot')) == 1)
    if max(strcmp(varargin, 'mean')) == 1 % Plot mean waveform
        plot(mUP);
    elseif max(strcmp(varargin, 'all')) == 1 % Plot all waveforms
        for kk = 1:length(E)
        plot(allUP(kk,:),'k');
        hold on;
        end
        hold on;
        plot(mUP,'r');
        hold off;
    end

    % Calculate features using the mean waveform
    pk = zeros(1,2);
    loc = zeros(1,2);
    pr = zeros(2,2);
    [pk2(1,1),pk2(2,1),w2,~] = findpeaks(mUP,'MinPeakHeight',20);
    [p,l,~,r] = findpeaks(mDOWN);

    if size(p,2) == 1 || size(l,2) == 1 || size(r,2) == 1
        pk(:,1) = p;
        loc(:,1) = l;
        pr(:,1)= r;
    elseif size(p,2) == 2
        pk = p;
        loc = l;
        pr = r;
    end

    F = mDOWN>0;
    I(1,1) = find(F,1,'first');
    I(1,2) = find(F,1,'last');
        % Params
        wf.a = pr(1,1);
        wf.b = pr(1,2);
        wf.pk1(1,1) = pk(1,1);
        wf.pk1(1,2) = loc(1,1);
        wf.pk2 = transpose(pk2);
        wf.pk3(1,1) = pk(1,2);
        wf.pk3(1,2) = loc(1,2);
        wf.c = wf.pk2(1,2) - wf.pk1(1,2);
        wf.d = wf.pk3(1,2) - wf.pk2(1,2);
        wf.e = wf.pk2(1,1) + wf.pk1(1,1);
        wf.f = wf.pk3(1,1) + wf.pk2(1,1);
        wf.g = w2;
        wf.h = wf.pk2(1,2) - I(1,1);
        wf.i = I(1,2) - wf.pk2(1,2);
        wf.mWF = mUP;

elseif max(strcmp(varargin, 'plot')) == 1
    if max(strcmp(varargin, 'mean')) == 1 % Plot mean waveform
        plot(mUP);
    elseif max(strcmp(varargin, 'all')) == 1 % Plot all waveforms
        for kk = 1:length(E)
        plot(allUP(kk,:),'k');
        hold on;
        end
        hold on;
        plot(mUP,'r');
        hold off;
    end
    
elseif max(strcmp(varargin, 'features')) == 1
    % Calculate features using the mean waveform
    pk = zeros(1,2);
    loc = zeros(1,2);
    pr = zeros(2,2);
    [pk2(1,1),pk2(2,1),w2,~] = findpeaks(mUP,'MinPeakHeight',20);
    [p,l,~,r] = findpeaks(mDOWN);

    if size(p,2) == 1 || size(l,2) == 1 || size(r,2) == 1
        pk(:,1) = p;
        loc(:,1) = l;
        pr(:,1)= r;
    elseif size(p,2) == 2
        pk = p;
        loc = l;
        pr = r;
    end

    F = mDOWN>0;
    I(1,1) = find(F,1,'first');
    I(1,2) = find(F,1,'last');
        % Params
        wf.a = pr(1,1);
        wf.b = pr(1,2);
        wf.pk1(1,1) = pk(1,1);
        wf.pk1(1,2) = loc(1,1);
        wf.pk2 = transpose(pk2);
        wf.pk3(1,1) = pk(1,2);
        wf.pk3(1,2) = loc(1,2);
        wf.c = wf.pk2(1,2) - wf.pk1(1,2);
        wf.d = wf.pk3(1,2) - wf.pk2(1,2);
        wf.e = wf.pk2(1,1) + wf.pk1(1,1);
        wf.f = wf.pk3(1,1) + wf.pk2(1,1);
        wf.g = w2;
        wf.h = wf.pk2(1,2) - I(1,1);
        wf.i = I(1,2) - wf.pk2(1,2);
        wf.mWF = mUP;
        
else % If there are no warargins plot mean waveform
    plot(mUP);
end

if s == 1
    saveas(gcf, genvarname([neuron(1:end-4),'_',num2str(G),'_wave']), type);
elseif s ~= 1
end    


end
