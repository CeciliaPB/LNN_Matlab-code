%Gergo Nagy, work in progress

[data, timestamps, info] = load_open_ephys_data('all_channels.events')

laserTTL = timestamps (1:6:end,:);
shockTTL = timestamps (5:6:end,:);
shockduringlaserTTL = timestamps (2:6:end,:);


allTTL = timestamps();