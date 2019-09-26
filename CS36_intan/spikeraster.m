function H = spikeraster(spiketimes,eventtimes,win)
%SPIKERASTER   Raster plot.
%   H = SPIKERASTER(SPIKETIMES,EVENTTIMES,WIN) plots one dot per spike
%   given in SPIKETIMES aligned to EVENTTIMES. A time window defined by WIN
%   is plotted. WIN should be a 1-by-2 array with a window start and end
%   relative to 0.

%   Balazs Hangya, TENSS 2016
%   hangya.balazs@koki.mta.hu

numEvents = length(eventtimes);

H = figure;
hold on
for iE = 1: numEvents
    lspk = spiketimes - eventtimes(iE);   % center spikes on the event
    lspk = lspk(lspk>win(1)&lspk<win(2));   % apply the window
    plot(lspk,iE*ones(size(lspk)),'k.','MarkerSize',12)   % plot
end