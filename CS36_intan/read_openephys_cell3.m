%%

DATAPATH = 'D:\Dropbox\courses\2015_TENSS\data\';
[data1, timestamps, info] = load_open_ephys_data([DATAPATH 'gezo_2015-06-07_18-12-18\100_CH5.continuous']);
[data2, timestamps, info] = load_open_ephys_data([DATAPATH 'gezo_2015-06-07_18-12-18\100_CH6.continuous']);
[data3, timestamps, info] = load_open_ephys_data([DATAPATH 'gezo_2015-06-07_18-12-18\100_CH7.continuous']);
[data4, timestamps, info] = load_open_ephys_data([DATAPATH 'gezo_2015-06-07_18-12-18\100_CH8.continuous']);
data = [data1 data2 data3 data4];

%%

[AllTimeStamps AllWaveForms] = intandisc3(data,info.header.sampleRate,70,[]);

%%

figure
plot(squeeze(AllWaveForms{1}(:,4,:))')


%%

avw = squeeze(mean(AllWaveForms{1},1));
figure;plot(avw')

%%

mx = squeeze(max(AllWaveForms{1},[],3));
mn = squeeze(min(AllWaveForms{1},[],3));
m = mx + mn;
figure;plot(m(:,1),mx(:,4),'.')