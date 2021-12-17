
function stats  = comp_PSTH_WNvsBSL_statszero(NumPartitions, spiketime_b,spiketime_t,statpsth_b,statpsth_t,dt,win_b,win_t,tags,sig_thr)
%comp_PSTH_WNvsBSL_statszero   Statistics for firing rate change.
%   The function tests whether spike number in a window of interest
%   change compared to a baseline period and, if it's the case, report it
%   in the raster plot produced.
%
%   See also comp_PSTH_WNvsBSL_stats and PSTH_STATS.
% -------------------------------------------------------------------------
% Modified from PSTH_STATS_WM_VS_BSLN by Nicola Solari
% 
% Cecília Pardo-Bellver, 2021
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
% -------------------------------------------------------------------------

% Default parameters
try
    [mylabels, mycolors] = makeColorsLabels(@defineLabelsColors_Cecilia,tags);
catch
    mylabels = tags;
    mycolors = [0 0.4470 0.7410];
end
testwin     = win_t;
baselinewin = win_b;
threshold   = 0.1;
ytx_multiplier1 = {0.8 0.6 0.4 0.2};
ytx_multiplier2 = {0.7 0.5 0.3 0.1};

if NumPartitions == 1
    % Specific spiketimes and psth
    spk_b  = spiketime_b;
    spk_t  = spiketime_t;
    psth_b = statpsth_b;
    psth_t = statpsth_t;
    
    % Index for time 0
    timew_b     = abs(win_b(1)) / dt;
    nullindex_b = timew_b + 1;
    timew_t     = abs(win_t(1)) / dt;
    nullindex_t = timew_t + 1;
    
    % Window for testing the potential effect
    WNb = [baselinewin(1)/dt+nullindex_b baselinewin(2)/dt+nullindex_b-1]; % baseline window; convert to indices; exlude 0
    WNt = [testwin(1)/dt+nullindex_t testwin(2)/dt+nullindex_t-1]; % test window; convert to indices
    WNb = round(WNb);
    WNt = round(WNt);
    lWNb = WNb(2) - WNb(1) + 1; % length of baseline window
    lWNt = WNt(2) - WNt(1) + 1; % length of test window
    % Trial number
    nTrial = size(spk_b,1);
    % Time vector
    time_b = win_b(1):dt:win_b(2);
    time_t = win_t(1):dt:win_t(2); % maybe?
    
    % Spiking probability -------------------------------------------------
    spkprob_b = spk_b(:,1:timew_b); % st < time 0  % ASSUMING BASELINE BEFORE TRIGGER ONE
    if testwin(1) < 0 % dealing with window position relatively to trigger event
        spkprob_t = spk_t(:,end - timew_t + 1  : end);
    else
        spkprob_t = spk_t(:,timew_t+1:end);
    end
    prob_b = sum(spkprob_b) / nTrial;   % spiking prob. (strictly) before time 0
    baseline_prob = mean(prob_b(WNb(1):WNb(2)))/dt;  % spikes/sec (was spikes/bin before)
    
    % Inhibition ----------------------------------------------------------
    minafter = min(psth_t(WNt(1):WNt(2)));
    if minafter > baseline_prob % unless it goes below baseline, no inhibition
        inhibition_start = NaN;
        inhibition_end   = NaN;
        inhibition_peak  = NaN;
        inhibition_time  = 0; % if firing does not go below baseline
        MWp_i = 0.99; % Mann-Whitney p-value, inhibition
    else % if it goes below baseline, check inhibition
        if testwin(1) < 0
            mininx = find(psth_t(WNt(1):WNt(2)) == minafter,1,'first'); % minimal firing
        else
            mininx = timew_t + find(psth_t(WNt(1):WNt(2)) == minafter,1,'first'); % minimal firing
        end
        thr       = baseline_prob - (baseline_prob - minafter) * threshold;  % threshold is determined in proportion of peak-baseline distance
        %         crossi_dw = valuecrossing(time_t(WNt(1):mininx),psth_t(WNt(1):mininx),thr,'down'); % Nicola's way
        %         crossi_dw_inx = valuecrossing(WNt(1):mininx,psth_t(WNt(1):mininx),thr,'down');
        timeline  = WNt(1):WNt(2)+1;
        dw = find(psth_t(WNt(1):mininx) >= thr & [psth_t(WNt(1)+1:mininx) thr] < thr);
        try
            crossi_dw_inx = interp1(psth_t(dw(1):dw(1)+1),(dw(1):dw(1)+1),thr);
            crossi_dw     = interp1(timeline,time_t,crossi_dw_inx);
        catch
            crossi_dw_inx = WNt(1);
            crossi_dw     = time_t(WNt(1));
        end
        if isempty(crossi_dw) || isnan(crossi_dw_inx)
            crossi_dw_inx = WNt(1);
            crossi_dw     = time_t(WNt(1));
        end
        crossi_dw_inx = round(crossi_dw_inx(end));
        %         crossi_up     = valuecrossing(time_t(mininx:WNt(2)),psth_t(mininx:WNt(2)),thr,'up'); % Nicola's way
        %         crossi_up_inx = valuecrossing(mininx:WNt(2),psth_t(mininx:WNt(2)),thr,'up');
        up = find(psth_t(mininx:WNt(2)) <= thr & [psth_t(mininx+1:WNt(2)) thr] > thr)+ mininx;
        try
            crossi_up_inx = interp1(psth_t(up(1)-1:up(1)),(up(1)-1:up(1)),thr);
            crossi_up     = interp1(timeline,time_t,crossi_up_inx);
        catch
            crossi_up_inx = WNt(2);
            crossi_up     = time_t(WNt(2));
        end
        if isempty(crossi_up) || isnan(crossi_up_inx)
            crossi_up_inx = WNt(2);
            crossi_up     = time_t(WNt(2));
        end
        crossi_up_inx    = round(crossi_up_inx(1));
        inhibition_start = crossi_dw(end); % last crossing of half-baseline probability before minimum
        inhibition_end   = crossi_up(1); % first crossing of half-baseline probability after minimum
        inhibition_time  = inhibition_end - inhibition_start;
        try
            [~,inhibition_peak]  = findpeaks((psth_t(crossi_dw_inx:crossi_up_inx))*-1,...
                time_t(crossi_dw_inx:crossi_up_inx),'NPeaks',1); % peak time of inhibition
        catch
            inhibition_peak = time_t(mininx) - time_t(timew_t+1);
        end
        % Nullhypothesis distribution
        ncross   = crossi_up_inx - crossi_dw_inx + 1;
        wlength  = floor(lWNb/lWNt); % split up the baseline according to test window length
        nulldist = nan(ncross,wlength*nTrial);
        for kk = 1:wlength
            inx = WNb(2)-kk*lWNt+1 : WNb(2)-(kk-1)*lWNt; % Baseline values, adjusted to test length
            spk_win  = spkprob_b(:,inx);
            psth_win = psth_b(inx);
            min_psth = find(psth_win == min(psth_win));
            min_psth = min_psth(1);
            inx2 = min_psth-floor(ncross/2) : min_psth+ceil(ncross/2)-1; % Baseline values, adjusted to ncross
            inx2 = inx2 - min(0,inx2(1)-1) - (max(length(inx),inx2(end)) - length(inx));
            if ~ismember(min_psth,inx2)
                error('vipisinfluenced2:nullhypoIndexing','Programming error.')
            end
            spk_win2 = spk_win(:,inx2);
            if ~isequal(size(spk_win2,2),ncross)
                error('vipisinfluenced2:nullhypoIndexing','Programming error.')
            end
            nulldist(:,(kk-1)*nTrial+1:kk*nTrial) = spk_win2';
        end
        if any(isnan(nulldist))
            error('vipisinfluenced2:nullhypoIndexing','Programming error.')
        end
        null_dist = sum(nulldist);
        % Test distribution
        if testwin(1) < 0
            test_dist = sum(spkprob_t(:,crossi_dw_inx : crossi_up_inx),2); % fuckoff, I give the correct window since the beginning, so half of these things are useless
        else
            test_dist = sum(spkprob_t(:,crossi_dw_inx-timew_t : crossi_up_inx-timew_t),2);
        end
        % Mann-Whitney test
        [MWp_i,MWh_i] = ranksum(null_dist,test_dist,'alpha',sig_thr,'tail','right'); % Inhibition x>y, 'right tail'
        ranks_i = tiedrank([null_dist test_dist']);
        tp      = length(null_dist);
        nullranksum_i = mean(ranks_i(1:tp));
        testranksum_i = mean(ranks_i(tp+1:end));
        if testranksum_i > nullranksum_i    % one-sided test
            MWp_i = 0.99;
            MWh_i = 0;
        end
        if MWh_i
            clri = mycolors; %[0 153 255] / 256;
            %     else
            %         clri =  [0 0 0];
        end
    end
    
    % Activation time -----------------------------------------------------
    maxafter = max(psth_t(WNt(1):WNt(2)));
    if maxafter < baseline_prob     % putative activation, if firing goes above baseline
        activation_start = NaN;
        activation_end   = NaN;
        activation_peak  = NaN;
        activation_time  = 0;   % if firing does not go above baseline
        MWp_a = 0.99;
    else
        if testwin(1) < 0
            maxinx = find(psth_t(WNt(1):WNt(2))==maxafter,1,'first');   % minimal firing
        else
            maxinx = timew_t + find(psth_t(WNt(1):WNt(2))==maxafter,1,'first');   % minimal firing
        end
        thr       = baseline_prob + (maxafter - baseline_prob) * threshold;  % threshold is determined in proportion of peak-baseline distance
        %         crossa_up = valuecrossing(time_t(WNt(1):maxinx),psth_t(WNt(1):maxinx),thr,'up');
        %         crossa_up_inx = valuecrossing(WNt(1):maxinx,psth_t(WNt(1):maxinx),thr,'up');
        timeline  = WNt(1):WNt(2)+1;
        up = find(psth_t(WNt(1):maxinx) <= thr & psth_t(WNt(1)+1:maxinx+1) > thr);
        try
            crossa_up_inx = interp1(psth_t(up(1):up(1)+1),(up(1):up(1)+1),thr);
            crossa_up     = interp1(timeline,time_t,crossi_up_inx);
        catch
            crossa_up_inx = WNt(1);
            crossa_up     = time_t(WNt(1));
        end
        if isempty(crossa_up) || isnan(crossa_up_inx)
            crossa_up_inx = WNt(1);
            crossa_up     = time_t(WNt(1));
        end
        crossa_up_inx = round(crossa_up_inx(end));
        %         crossa_dw     = valuecrossing(time_t(maxinx:WNt(2)),psth_t(maxinx:WNt(2)),thr,'down');
        %         crossa_dw_inx = valuecrossing(maxinx:WNt(2),psth_t(maxinx:WNt(2)),thr,'down');
        dw = find(psth_t(maxinx:WNt(2)) >= thr & [psth_t(maxinx+1:WNt(2)) thr] <= thr)+maxinx;
        try
            crossa_dw_inx = interp1(psth_t(dw(1)-1:dw(1)),(dw(1)-1:dw(1)),thr);
            crossa_dw     = interp1(timeline,time_t,crossi_dw_inx);
        catch
            crossa_dw_inx = WNt(2);
            crossa_dw     = time_t(WNt(2));
        end
        if isempty(crossa_dw) || isnan(crossa_dw_inx)
            crossa_dw_inx = WNt(2);
            crossa_dw     = time_t(WNt(2));
        end
        crossa_dw_inx    = round(crossa_dw_inx(1));
        activation_start = crossa_up(end);   % last crossing of one and a half-baseline probability before maximum
        activation_end   = crossa_dw(1);   % first crossing of one and a half-baseline probability after maximum
        activation_time  = activation_end - activation_start;
        try
            [~,activation_peak]  = findpeaks((psth_t(crossa_up_inx:crossa_dw_inx)),...
                time_t(crossa_up_inx:crossa_dw_inx),'NPeaks',1,'SortStr','descend'); % peak time of inhibition
        catch
            activation_peak = time_t(maxinx) - time_t(timew_t+1);
        end
        % Nullhypothesis distribution
        ncross   = crossa_dw_inx - crossa_up_inx + 1;
        wlength  = floor(lWNb/lWNt);   % split up the baseline according to test window length
        nulldist = nan(ncross,wlength*nTrial);
        for kk = 1:wlength
            inx = WNb(2)-kk*lWNt+1:WNb(2)-(kk-1)*lWNt;
            spk_win  = spkprob_b(:,inx);
            psth_win = psth_b(inx);
            min_psth = find(psth_win==max(psth_win));
            min_psth = min_psth(1);
            inx2 = min_psth - floor(ncross/2) : min_psth + ceil(ncross/2) - 1;
            inx2 = inx2 - min(0,inx2(1)-1) - (max(length(inx),inx2(end)) - length(inx));
            if ~ismember(min_psth,inx2)
                error('vipisinfluenced2:nullhypoIndexing','Programming error.')
            end
            spk_win2 = spk_win(:,inx2);
            if ~isequal(size(spk_win2,2),ncross)
                error('vipisinfluenced2:nullhypoIndexing','Programming error.')
            end
            nulldist(:,(kk-1)*nTrial+1:kk*nTrial) = spk_win2';
        end
        if any(isnan(nulldist))
            error('vipisinfluenced2:nullhypoIndexing','Programming error.')
        end
        null_dist = sum(nulldist);
        % Test distribution
        if testwin(1) < 0
            test_dist = sum(spkprob_t(:,crossa_up_inx:crossa_dw_inx),2);
        else
            test_dist = sum(spkprob_t(:,crossa_up_inx-timew_t:crossa_dw_inx-timew_t),2);
        end
        % Mann-Whitney test
        [MWp_a,MWh_a] = ranksum(null_dist,test_dist,'alpha',sig_thr,'tail','left'); % Activation x<y, 'left tail'
        ranks_a = tiedrank([null_dist test_dist']);
        tp      = length(null_dist);
        nullranksum_a = mean(ranks_a(1:tp));
        testranksum_a = mean(ranks_a(tp+1:end));
        if testranksum_a < nullranksum_a    % one-sided test
            MWp_a = 0.99;
            MWh_a = 0;
        end
        if MWh_a
            clra = mycolors; % 'red'
        end
    end
    
    % % Plot
    if testwin(1) < 0
        h = plot([time_b abs(flip(time_t))],[psth_b psth_t],'color',mycolors);
        xlabel('Time (s)');
        ylabel('Rate (Hz)');
        xlim([time_b(1) abs(time_t(1))])
    else
        h =plot([time_b time_t],[psth_b psth_t],'color',mycolors);
        %         legend(mylabels,'Box','off','Location','northwest');
        xlabel('Time (s)');
        ylabel('Rate (Hz)');
        xlim([time_b(1) time_t(end)]);
    end
    
    if exist('clri','var')
        hold on
        if testwin(1) < 0
            newt   = abs(flip(time_t));
            hI = plot(newt(crossi_dw_inx:crossi_up_inx),psth_t(crossi_dw_inx:crossi_up_inx),'color',...
                mycolors,'LineWidth',2);
        else
            hI = plot(time_t(crossi_dw_inx:crossi_up_inx),psth_t(crossi_dw_inx:crossi_up_inx),'color',...
                mycolors,'LineWidth',2);
        end
        x_lim = xlim;
        y_lim = ylim;
        text(x_lim(1)+(x_lim(2)-x_lim(1))*0.4,y_lim(1)+(y_lim(2)-y_lim(1))*ytx_multiplier1{1},...
            ['\itInhib p = ',num2str(MWp_i)],'Color',mycolors);
    end
    if exist('clra','var')
        hold on
        if testwin(1) < 0
            newt   = abs(flip(time_t));
            hA = plot(newt(crossa_up_inx:crossa_dw_inx),psth_t(crossa_up_inx:crossa_dw_inx),'color',...
                mycolors,'LineWidth',2);
        else
            hA = plot(time_t(crossa_up_inx:crossa_dw_inx),psth_t(crossa_up_inx:crossa_dw_inx),'color',...
                mycolors,'LineWidth',2);
        end
        x_lim = xlim;
        y_lim = ylim;
        text(x_lim(1)+(x_lim(2)-x_lim(1))*0.4,y_lim(1)+(y_lim(2)-y_lim(1))*ytx_multiplier2{1},...
            ['\itExcit  p = ',num2str(MWp_a)],'Color',mycolors);
    end
    hold on
    
    % Output statistics
    stats.mean_baseline     = mean(psth_b(WNb(1):WNb(2))); %#ok<*AGROW>
    stats.mean_test         = mean(psth_t(WNt(1):WNt(2)));
    stats.baseline          = baseline_prob;
    stats.minvalue          = minafter;
    stats.inhibition_start  = inhibition_start;
    stats.inhibition_end    = inhibition_end;
    stats.inhibition_peak   = inhibition_peak;
    stats.inhibition_time   = inhibition_time;
    stats.MWp_i             = MWp_i;
    stats.maxvalue          = maxafter;
    stats.activation_start  = activation_start;
    stats.activation_end    = activation_end;
    stats.activation_peak   = activation_peak;
    stats.activation_time   = activation_time;
    stats.MWp_a             = MWp_a;
    stats.tag               = mylabels;
    
    legend(h,mylabels,'Box','off','Location','northwest');
    
elseif NumPartitions > 1
    for iP = 1:NumPartitions
        % Specific spiketimes and psth
        spk_b  = spiketime_b{iP};
        spk_t  = spiketime_t{iP};
        psth_b = statpsth_b(iP,:);
        psth_t = statpsth_t(iP,:);
        
        % Index for time 0
        timew_b     = abs(win_b(1)) / dt;
        nullindex_b = timew_b + 1;
        timew_t     = abs(win_t(1)) / dt;
        nullindex_t = timew_t + 1;
        
        % Window for testing the potential effect
        WNb = [baselinewin(1)/dt+nullindex_b baselinewin(2)/dt+nullindex_b-1]; % baseline window; convert to indices; exlude 0
        WNt = [testwin(1)/dt+nullindex_t testwin(2)/dt+nullindex_t-1]; % test window; convert to indices
        WNb = round(WNb);
        WNt = round(WNt);
        lWNb = WNb(2) - WNb(1) + 1; % length of baseline window
        lWNt = WNt(2) - WNt(1) + 1; % length of test window
        % Trial number
        nTrial = size(spk_b,1);
        % Time vector
        time_b = win_b(1):dt:win_b(2);
        time_t = win_t(1):dt:win_t(2); % maybe?
        
        % Spiking probability -------------------------------------------------
        spkprob_b = spk_b(:,1:timew_b); % st < time 0  % ASSUMING BASELINE BEFORE TRIGGER ONE
        if testwin(1) < 0 % dealing with window position relatively to trigger event
            spkprob_t = spk_t(:,end - timew_t + 1  : end);
        else
            spkprob_t = spk_t(:,timew_t+1:end);
        end
        prob_b = sum(spkprob_b) / nTrial;   % spiking prob. (strictly) before time 0
        baseline_prob = mean(prob_b(WNb(1):WNb(2)))/dt;  % spikes/sec (was spikes/bin before)
        
        % Inhibition ----------------------------------------------------------
        minafter = min(psth_t(WNt(1):WNt(2)));
        if minafter > baseline_prob % unless it goes below baseline, no inhibition
            inhibition_start = NaN;
            inhibition_end   = NaN;
            inhibition_peak  = NaN;
            inhibition_time  = 0; % if firing does not go below baseline
            MWp_i = 0.99; % Mann-Whitney p-value, inhibition
        else % if it goes below baseline, check inhibition
            if testwin(1) < 0
                mininx = find(psth_t(WNt(1):WNt(2)) == minafter,1,'first'); % minimal firing
            else
                mininx = timew_t + find(psth_t(WNt(1):WNt(2)) == minafter,1,'first'); % minimal firing
            end
            thr       = baseline_prob - (baseline_prob - minafter) * threshold;  % threshold is determined in proportion of peak-baseline distance
            %         crossi_dw = valuecrossing(time_t(WNt(1):mininx),psth_t(WNt(1):mininx),thr,'down'); % Nicola's way
            %         crossi_dw_inx = valuecrossing(WNt(1):mininx,psth_t(WNt(1):mininx),thr,'down');
            timeline  = WNt(1):WNt(2)+1;
            dw = find(psth_t(WNt(1):mininx) >= thr & [psth_t(WNt(1)+1:mininx) thr] < thr);
            try
                crossi_dw_inx = interp1(psth_t(dw(1):dw(1)+1),(dw(1):dw(1)+1),thr);
                crossi_dw     = interp1(timeline,time_t,crossi_dw_inx);
            catch
                crossi_dw_inx = WNt(1);
                crossi_dw     = time_t(WNt(1));
            end
            if isempty(crossi_dw) || isnan(crossi_dw_inx)
                crossi_dw_inx = WNt(1);
                crossi_dw     = time_t(WNt(1));
            end
            crossi_dw_inx = round(crossi_dw_inx(end));
            %         crossi_up     = valuecrossing(time_t(mininx:WNt(2)),psth_t(mininx:WNt(2)),thr,'up'); % Nicola's way
            %         crossi_up_inx = valuecrossing(mininx:WNt(2),psth_t(mininx:WNt(2)),thr,'up');
            up = find(psth_t(mininx:WNt(2)) <= thr & [psth_t(mininx+1:WNt(2)) thr] > thr)+ mininx;
            try
                crossi_up_inx = interp1(psth_t(up(1)-1:up(1)),(up(1)-1:up(1)),thr);
                crossi_up     = interp1(timeline,time_t,crossi_up_inx);
            catch
                crossi_up_inx = WNt(2);
                crossi_up     = time_t(WNt(2));
            end
            if isempty(crossi_up) || isnan(crossi_up_inx)
                crossi_up_inx = WNt(2);
                crossi_up     = time_t(WNt(2));
            end
            crossi_up_inx    = round(crossi_up_inx(1));
            inhibition_start = crossi_dw(end); % last crossing of half-baseline probability before minimum
            inhibition_end   = crossi_up(1); % first crossing of half-baseline probability after minimum
            inhibition_time  = inhibition_end - inhibition_start;
            try
                [~,inhibition_peak]  = findpeaks((psth_t(crossi_dw_inx:crossi_up_inx))*-1,...
                    time_t(crossi_dw_inx:crossi_up_inx),'NPeaks',1); % peak time of inhibition
            catch
                inhibition_peak = time_t(mininx) - time_t(timew_t+1);
            end
            % Nullhypothesis distribution
            ncross   = crossi_up_inx - crossi_dw_inx + 1;
            wlength  = floor(lWNb/lWNt); % split up the baseline according to test window length
            nulldist = nan(ncross,wlength*nTrial);
            for kk = 1:wlength
                inx = WNb(2)-kk*lWNt+1 : WNb(2)-(kk-1)*lWNt; % Baseline values, adjusted to test length
                spk_win  = spkprob_b(:,inx);
                psth_win = psth_b(inx);
                min_psth = find(psth_win == min(psth_win));
                min_psth = min_psth(1);
                inx2 = min_psth-floor(ncross/2) : min_psth+ceil(ncross/2)-1; % Baseline values, adjusted to ncross
                inx2 = inx2 - min(0,inx2(1)-1) - (max(length(inx),inx2(end)) - length(inx));
                if ~ismember(min_psth,inx2)
                    error('vipisinfluenced2:nullhypoIndexing','Programming error.')
                end
                spk_win2 = spk_win(:,inx2);
                if ~isequal(size(spk_win2,2),ncross)
                    error('vipisinfluenced2:nullhypoIndexing','Programming error.')
                end
                nulldist(:,(kk-1)*nTrial+1:kk*nTrial) = spk_win2';
            end
            if any(isnan(nulldist))
                error('vipisinfluenced2:nullhypoIndexing','Programming error.')
            end
            null_dist = sum(nulldist);
            % Test distribution
            if testwin(1) < 0
                test_dist = sum(spkprob_t(:,crossi_dw_inx : crossi_up_inx),2); % fuckoff, I give the correct window since the beginning, so half of these things are useless
            else
                test_dist = sum(spkprob_t(:,crossi_dw_inx-timew_t : crossi_up_inx-timew_t),2);
            end
            % Mann-Whitney test
            [MWp_i,MWh_i] = ranksum(null_dist,test_dist,'alpha',sig_thr,'tail','right'); % Inhibition x>y, 'right tail'
            ranks_i = tiedrank([null_dist test_dist']);
            tp      = length(null_dist);
            nullranksum_i = mean(ranks_i(1:tp));
            testranksum_i = mean(ranks_i(tp+1:end));
            if testranksum_i > nullranksum_i    % one-sided test
                MWp_i = 0.99;
                MWh_i = 0;
            end
            if MWh_i
                clri = mycolors{iP}; %[0 153 255] / 256;
                %     else
                %         clri =  [0 0 0];
            end
        end
        
        % Activation time -----------------------------------------------------
        maxafter = max(psth_t(WNt(1):WNt(2)));
        if maxafter < baseline_prob     % putative activation, if firing goes above baseline
            activation_start = NaN;
            activation_end   = NaN;
            activation_peak  = NaN;
            activation_time  = 0;   % if firing does not go above baseline
            MWp_a = 0.99;
        else
            if testwin(1) < 0
                maxinx = find(psth_t(WNt(1):WNt(2))==maxafter,1,'first');   % minimal firing
            else
                maxinx = timew_t + find(psth_t(WNt(1):WNt(2))==maxafter,1,'first');   % minimal firing
            end
            thr       = baseline_prob + (maxafter - baseline_prob) * threshold;  % threshold is determined in proportion of peak-baseline distance
            %         crossa_up = valuecrossing(time_t(WNt(1):maxinx),psth_t(WNt(1):maxinx),thr,'up');
            %         crossa_up_inx = valuecrossing(WNt(1):maxinx,psth_t(WNt(1):maxinx),thr,'up');
            timeline  = WNt(1):WNt(2)+1;
            up = find(psth_t(WNt(1):maxinx) <= thr & psth_t(WNt(1)+1:maxinx+1) > thr);
            try
                crossa_up_inx = interp1(psth_t(up(1):up(1)+1),(up(1):up(1)+1),thr);
                crossa_up     = interp1(timeline,time_t,crossi_up_inx);
            catch
                crossa_up_inx = WNt(1);
                crossa_up     = time_t(WNt(1));
            end
            if isempty(crossa_up) || isnan(crossa_up_inx)
                crossa_up_inx = WNt(1);
                crossa_up     = time_t(WNt(1));
            end
            crossa_up_inx = round(crossa_up_inx(end));
            %         crossa_dw     = valuecrossing(time_t(maxinx:WNt(2)),psth_t(maxinx:WNt(2)),thr,'down');
            %         crossa_dw_inx = valuecrossing(maxinx:WNt(2),psth_t(maxinx:WNt(2)),thr,'down');
            dw = find(psth_t(maxinx:WNt(2)) >= thr & [psth_t(maxinx+1:WNt(2)) thr] <= thr)+maxinx;
            try
                crossa_dw_inx = interp1(psth_t(dw(1)-1:dw(1)),(dw(1)-1:dw(1)),thr);
                crossa_dw     = interp1(timeline,time_t,crossi_dw_inx);
            catch
                crossa_dw_inx = WNt(2);
                crossa_dw     = time_t(WNt(2));
            end
            if isempty(crossa_dw) || isnan(crossa_dw_inx)
                crossa_dw_inx = WNt(2);
                crossa_dw     = time_t(WNt(2));
            end
            crossa_dw_inx    = round(crossa_dw_inx(1));
            activation_start = crossa_up(end);   % last crossing of one and a half-baseline probability before maximum
            activation_end   = crossa_dw(1);   % first crossing of one and a half-baseline probability after maximum
            activation_time  = activation_end - activation_start;
            try
                [~,activation_peak]  = findpeaks((psth_t(crossa_up_inx:crossa_dw_inx)),...
                    time_t(crossa_up_inx:crossa_dw_inx),'NPeaks',1,'SortStr','descend'); % peak time of inhibition
            catch
                activation_peak = time_t(maxinx) - time_t(timew_t+1);
            end
            % Nullhypothesis distribution
            ncross   = crossa_dw_inx - crossa_up_inx + 1;
            wlength  = floor(lWNb/lWNt);   % split up the baseline according to test window length
            nulldist = nan(ncross,wlength*nTrial);
            for kk = 1:wlength
                inx = WNb(2)-kk*lWNt+1:WNb(2)-(kk-1)*lWNt;
                spk_win  = spkprob_b(:,inx);
                psth_win = psth_b(inx);
                min_psth = find(psth_win==max(psth_win));
                min_psth = min_psth(1);
                inx2 = min_psth - floor(ncross/2) : min_psth + ceil(ncross/2) - 1;
                inx2 = inx2 - min(0,inx2(1)-1) - (max(length(inx),inx2(end)) - length(inx));
                if ~ismember(min_psth,inx2)
                    error('vipisinfluenced2:nullhypoIndexing','Programming error.')
                end
                spk_win2 = spk_win(:,inx2);
                if ~isequal(size(spk_win2,2),ncross)
                    error('vipisinfluenced2:nullhypoIndexing','Programming error.')
                end
                nulldist(:,(kk-1)*nTrial+1:kk*nTrial) = spk_win2';
            end
            if any(isnan(nulldist))
                error('vipisinfluenced2:nullhypoIndexing','Programming error.')
            end
            null_dist = sum(nulldist);
            % Test distribution
            if testwin(1) < 0
                test_dist = sum(spkprob_t(:,crossa_up_inx:crossa_dw_inx),2);
            else
                test_dist = sum(spkprob_t(:,crossa_up_inx-timew_t:crossa_dw_inx-timew_t),2);
            end
            % Mann-Whitney test
            [MWp_a,MWh_a] = ranksum(null_dist,test_dist,'alpha',sig_thr,'tail','left'); % Activation x<y, 'left tail'
            ranks_a = tiedrank([null_dist test_dist']);
            tp      = length(null_dist);
            nullranksum_a = mean(ranks_a(1:tp));
            testranksum_a = mean(ranks_a(tp+1:end));
            if testranksum_a < nullranksum_a    % one-sided test
                MWp_a = 0.99;
                MWh_a = 0;
            end
            if MWh_a
                clra = mycolors{iP}; % 'red'
            end
        end
        
        % % Plot
        if testwin(1) < 0
            h(iP) = plot([time_b abs(flip(time_t))],[psth_b psth_t],'color',mycolors{iP});
            xlabel('Time (s)');
            ylabel('Rate (Hz)');
            xlim([time_b(1) abs(time_t(1))])
        else
            h(iP) =plot([time_b time_t],[psth_b psth_t],'color',mycolors{iP});
            %         legend(mylabels,'Box','off','Location','northwest');
            xlabel('Time (s)');
            ylabel('Rate (Hz)');
            xlim([time_b(1) time_t(end)]);
        end
        
        if exist('clri','var')
            hold on
            if testwin(1) < 0
                newt   = abs(flip(time_t));
                hI(iP) = plot(newt(crossi_dw_inx:crossi_up_inx),...
                    psth_t(crossi_dw_inx:crossi_up_inx),'color',...
                    mycolors{iP},'LineWidth',2);
            else
                hI(iP) = plot(time_t(crossi_dw_inx:crossi_up_inx),...
                    psth_t(crossi_dw_inx:crossi_up_inx),'color',...
                    mycolors{iP},'LineWidth',2);
            end
            x_lim = xlim;
            y_lim = ylim;
            text(x_lim(1)+(x_lim(2)-x_lim(1))*0.4,...
                y_lim(1)+(y_lim(2)-y_lim(1))*ytx_multiplier1{iP},...
                ['\itInhib p = ',num2str(MWp_i)],'Color',mycolors{iP});
        end
        if exist('clra','var')
            hold on
            if testwin(1) < 0
                newt   = abs(flip(time_t));
                hA(iP) = plot(newt(crossa_up_inx:crossa_dw_inx),...
                    psth_t(crossa_up_inx:crossa_dw_inx),'color',...
                    mycolors{iP},'LineWidth',2);
            else
                hA(iP) = plot(time_t(crossa_up_inx:crossa_dw_inx),...
                    psth_t(crossa_up_inx:crossa_dw_inx),'color',...
                    mycolors{iP},'LineWidth',2);
            end
            x_lim = xlim;
            y_lim = ylim;
            text(x_lim(1)+(x_lim(2)-x_lim(1))*0.4,...
                y_lim(1)+(y_lim(2)-y_lim(1))*ytx_multiplier2{iP},...
                ['\itExcit  p = ',num2str(MWp_a)],'Color',mycolors{iP});
        end
        hold on
        
        % Output statistics
        stats{iP}.mean_baseline     = mean(psth_b(WNb(1):WNb(2))); %#ok<*AGROW>
        stats{iP}.mean_test         = mean(psth_t(WNt(1):WNt(2)));
        stats{iP}.baseline          = baseline_prob;
        stats{iP}.minvalue          = minafter;
        stats{iP}.inhibition_start  = inhibition_start;
        stats{iP}.inhibition_end    = inhibition_end;
        stats{iP}.inhibition_peak   = inhibition_peak;
        stats{iP}.inhibition_time   = inhibition_time;
        stats{iP}.MWp_i             = MWp_i;
        stats{iP}.maxvalue          = maxafter;
        stats{iP}.activation_start  = activation_start;
        stats{iP}.activation_end    = activation_end;
        stats{iP}.activation_peak   = activation_peak;
        stats{iP}.activation_time   = activation_time;
        stats{iP}.MWp_a             = MWp_a;
        
    end
    legend(h,mylabels,'Box','off','Location','northwest');
end
