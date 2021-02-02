function [psth_spx, psth_t] = psth_hist(psth, bins)

% Generates a matrix for raster plot from cell array 'trialspx' generated
% by main function: ttl_psth.
%
% INPUTS 
% - psth: peri-stimulus time histogram
% - bins: kind like the resolution. 
% 
% OUTPUT
% - psth_spx: bars of the psth
% - psth_t: time of the psth
% 
% Examples 
% [psth_spx, psth_t] = psth_hist(psth, 1000)
%
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2019
% Laboratory of Network Neurophysiology
% Institute of Experimantal Medicine, Hungary.
% -------------------------------------------------------------------------

psth_spx = [];
psth_t = transpose(min(psth(:,1)) : bins : max(psth(:,1)));
spikes = psth(:,2);

  for ii = 1:bins:ceil(length(psth(:,1)))
      try 
          if ii < ceil(length(spikes))
            j = ii+bins-1;
            A = sum(spikes(ii:j));
            psth_spx = [psth_spx; A];
          elseif ii == ceil(length(spikes(:,1)))
            psth_spx = [psth_spx;spikes(ii)];
          end
      catch
          psth_spx = [psth_spx;spikes(ii)];
      end
  end 
end