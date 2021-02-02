
ChID = '101_CH20'; % Raw trace
neuID = [4,5]; % Neurons clustered in this trace [1,5:6]

% PARAMS ------------------------------------------------------------------
fr = 30000; % sampling frequency rate

% Paths 
oldpath = pwd;
neuronpath = [oldpath filesep 'GR'];

% Select group number
CHnum = str2double(ChID(7:end));
if CHnum <= 4
    GRnum = '1'; 
elseif CHnum > 4  && CHnum <= 8
    GRnum = '2';
elseif CHnum > 8  && CHnum <= 12
    GRnum = '3';
elseif CHnum > 12 && CHnum <= 16
    GRnum = '4';
elseif CHnum > 16 && CHnum <= 20
    GRnum = '5';
elseif CHnum > 20 && CHnum <= 24
    GRnum = '6';
elseif CHnum > 24 && CHnum <= 28
    GRnum = '7';
elseif CHnum > 28 && CHnum <= 32
    GRnum = '8';
end

% ColorCode up to 12 neurons 
%       TAG      LabelName              Color               
colorsID = { ...
        '1',     ['GR' GRnum ' - 1'],     [0.9 0.1 0];    ...
        '2',     ['GR' GRnum ' - 2'],     [0.9 0.6 0];    ...
        '3',     ['GR' GRnum ' - 3'],     [0.9 0.9 0.3];  ...
        '4',     ['GR' GRnum ' - 4'],     [0.9 0.3 0.4];  ...
        '5',     ['GR' GRnum ' - 5'],     [0.3 0.8 1];    ...
        '6',     ['GR' GRnum ' - 6'],     [0.4 0   1];    ...
        '7',     ['GR' GRnum ' - 7'],     [0.4 1   0.2];  ...
        '8',     ['GR' GRnum ' - 8'],     [0.2 0.6 0];    ...
        '9',     ['GR' GRnum ' - 9'],     [0.4 1   1];    ...
        '10',    ['GR' GRnum ' - 10'],    [0.3 0.6 1];    ...
        '11',    ['GR' GRnum ' - 11'],    [0.3 0.3 1];    ...
        '12',    ['GR' GRnum ' - 12'],    [0.8 0.1 0];    ...
        };
% -------------------------------------------------------------------------
 
% Load and plot channel data 
if exist('common_avg','var')== 1
elseif exist('common_avg','var')== 0
common_avg = common_avg_ref_probe(oldpath,32,[],'processor',str2double(ChID(1:3)),'rawdatafiletag','');
end
[data, ts] = load_open_ephys_data([ChID '.continuous']);
data = data - common_avg(:,1);
[b,a] = butter(3,[700 7000]/(fr/2),'bandpass');
fdata = filter(b,a,data);
% Plot filtered trace
subplot(2,1,2);
plot(ts,fdata,'k');
xlim([min(ts) max(ts)]);

% Change current directory and load neurons
cd(neuronpath);

for ii = neuID
    load(['GR' GRnum '_' num2str(ii) '.mat'],'TS');

    TT = ((TS(:,1)/10000));
    % Waveform window
    win = [-0.0003 0.0006];   % -300 to 600 us, same as used by MClust
    winp = round(win*fr);

    ts2 = nan(length(TT)*28,1);
    data2 = nan(length(TT)*28,1);
    color = colorsID{ii,3};
    A = find(neuID == ii);
    label(A,:) = colorsID{ii,2};
    for kk = 1:length(TT)
        timeSPK = find(ts == TT(kk)); 
        timeWF = [timeSPK+winp(1) timeSPK+winp(2)];
        jj = 27*kk + kk;
        ts2(jj-27:jj,1) = ts(timeWF(1):timeWF(2),1);
        data2(jj-27:jj,1) = fdata(timeWF(1):timeWF(2),1);   
    end
    
    subplot(2,1,2);
    hold on;
    for kk = 1:length(TT)
        jj = 27*kk + kk;
        plot(ts2(jj-27:jj,1),data2(jj-27:jj,1),'Color',color,'LineWidth',2);   
    end
    hold off;

    subplot(2,1,1);
    hold on;
    for kk = 1:length(TT)
        jj = 27*kk + kk;
        plot(ts2(jj-27:jj,1),data2(jj-27:jj,1),'Color',color,'LineWidth',2);   
    end
    legend(label);
    xlim([min(ts) max(ts)]);
    hold off;

end

cd(oldpath);
% print(gcf,'-depsc','-painters','figure.eps'); % save with vectorial
% % info