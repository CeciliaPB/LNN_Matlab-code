% plot_rawdata_wf To plot the raw data (common average referenced and
% filtered) and over it the selected waveforms. 
%
% IN PROGRESS
%
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2019
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
% -------------------------------------------------------------------------

[ChID, oldpath] = uigetfile('*.continuous*', 'Select Raw data file', 'MultiSelect', 'off'); % Raw trace
[neuID, neuronpath] = uigetfile('*.mat*', 'Select Neuron(s) data file', 'MultiSelect', 'on'); % Neurons timestamps

% PARAMS ------------------------------------------------------------------
fr = 30000; % sampling frequency rate

% Select group number
CHnum = str2double(ChID(7:end));
GRnum = floor(CHnum/4)+1;

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
 
% Load and plot channel data ----------------------------------------------
cd(oldpath)

common_avg = common_avg_ref_probe(oldpath,32,[],'processor',str2double(ChID(1:3)),...
    'rawdatafiletag','');
[data, ts] = load_open_ephys_data([ChID '.continuous']);
data = data - common_avg(:,1);
[b,a] = butter(3,[700 7000]/(fr/2),'bandpass');
fdata = filter(b,a,data);

% Plot filtered trace
subplot(2,1,2);
plot(ts,fdata,'k');
xlim([min(ts) max(ts)]);

% Change current directory and load neurons -------------------------------
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

% print(gcf,'-depsc','-painters','figure.eps'); % save with vectorial
% % info