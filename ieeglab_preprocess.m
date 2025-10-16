function EEG = ieeglab_preprocess(EEG)

% GUI to get user choices
if nargin < 2
    [EEG, wasCancelled] = ieeglab_gui_preprocess(EEG);
    if wasCancelled
        return  % user aborted, exit gracefully
    end
end

opt = EEG.ieeglab.opt;

% Keep only events of interest (print how many are removed)
if isfield(opt, 'event_filters')
    ev_tbl     = opt.events;
    ev_choices = opt.event_filters;
    vars       = fields(ev_choices);

    for iField = 1:numel(vars)
        labelToKeep = ev_choices.(vars{iField});
        if isempty(labelToKeep), continue; end

        % Remove rows whose value in this column is NOT in labelToKeep
        idxToRem = ~ismissing(ev_tbl.(vars{iField}), labelToKeep);

        if any(idxToRem)
            % Build readable list for printing
            vals = string(labelToKeep(:));
            vals = vals(~ismissing(vals));                 % drop <missing>
            labelStr = strjoin(cellstr(vals), ', ');       % "A, B, C"

            fprintf('Removing %d events that are not %s: {%s}\n', ...
                sum(idxToRem), char(vars{iField}), labelStr);

            % Apply removal to both table and EEG.event
            ev_tbl(idxToRem,:) = [];
            EEG.event(idxToRem) = [];
        end
    end
end
EEG.ieeglab.opt = rmfield(EEG.ieeglab.opt, 'events');

% Remove electrodes with no coordinates (if selected; print how many are removed)
if isfield(opt, 'remove_no_coords') && opt.remove_no_coords
    
    % keep only channels with valid XYZ
    hasXYZ = arrayfun(@(c) isfield(c,'X') && isfield(c,'Y') && isfield(c,'Z') && ...
        ~isempty(c.X) && ~isempty(c.Y) && ~isempty(c.Z) && ...
        all(isfinite([c.X c.Y c.Z])), EEG.chanlocs);
    keepChanIdx = find(hasXYZ);
    removedLabels = {EEG.chanlocs(~hasXYZ).labels};  % for logging
    if ~isempty(removedLabels)
        warning('Removing %d channels with no coordinates: %s\n', ...
            numel(removedLabels), strjoin(removedLabels, ', '));
    end
    EEG = pop_select(EEG, 'channel', keepChanIdx);
    EEG = eeg_checkset(EEG);

    % Remove corresponding events
    remEv = contains(lower({EEG.event.type}), lower(removedLabels));
    if any(remEv)
        warning("Removing %g events containing these electrodes with no coordinates.", sum(remEv))
        EEG.event(remEv) = [];
        EEG = eeg_checkset(EEG);
    end
end

% Remove events with electrodes not present in EEG.chanlocs
idx = ~contains(lower({EEG.event.type}), lower({EEG.chanlocs.labels}));
warning("Removing %g events containing an electrode that is not present in EEG.chanlocs: ", sum(idx))
disp(unique({EEG.event(idx).type}))
EEG.event(idx) = [];

% Downsample
if isfield(opt, 'downsample') && ~isempty(opt.downsample) && opt.downsample<EEG.srate
    fprintf("Downsampling iEEG data to %g Hz... \n", opt.downsample)
    EEG = pop_resample(EEG, opt.downsample);
end

% High-pass filter
if isfield(opt,'apply_highpass') && opt.apply_highpass && isfield(opt,'highpass') ...
        && ~isempty(opt.highpass) && opt.highpass > 0
    if isfield(opt,'highpass_type') && opt.highpass_type == 2
        EEG = pop_eegfiltnew(EEG, 'locutoff', double(opt.highpass), 'usefftfilt', 1, 'minphase', 1);
    else
        EEG = pop_eegfiltnew(EEG, 'locutoff', double(opt.highpass), 'usefftfilt', 1);
    end
end

% Notch filter
% Supports scalar or vector centers in choices.notch (e.g., 60 or [60 120 180])
if isfield(opt,'apply_notch') && opt.apply_notch && isfield(opt,'notch') && ~isempty(opt.notch)
    centers = double(opt.notch(:))';
    BW = 2;   % default ±1 Hz around each center; tweak here if you prefer another BW
    for f0 = centers
        lo = max(0, f0 - BW/2);
        hi = f0 + BW/2;
        EEG = pop_eegfiltnew(EEG, 'locutoff', lo, 'hicutoff', hi, 'usefftfilt', 1, 'revfilt', 1);
    end
end

% Lowpass filter
if isfield(opt,'apply_lowpass') && opt.apply_lowpass && isfield(opt,'lowpass') ...
        && ~isempty(opt.lowpass) && opt.lowpass > 0
    EEG = pop_eegfiltnew(EEG, 'hicutoff', double(opt.lowpass), 'usefftfilt', 1);
end

% Remove conditions with too few trials
if isfield(opt,'remove_rare_cond') && opt.remove_rare_cond && isfield(opt,'min_trials')
    minN = max(0, round(double(opt.min_trials)));
    if minN>0 && isfield(EEG,'event') && ~isempty(EEG.event)
        types = string({EEG.event.type});
        u = unique(types);
        cnt = arrayfun(@(x) sum(types==x), u);
        rm = u(cnt < minN);
        if ~isempty(rm)
            fprintf('Removing %d condition(s) with < %d trials: %s\n', numel(rm), minN, strjoin(cellstr(rm), ', '));
            keep = ~ismember(types, rm);
            EEG.event = EEG.event(keep);
            try EEG = eeg_checkset(EEG,'makeur'); catch, end
        else
            fprintf('No conditions below %d trials.\n', minN);
        end
    end
end


% Epoch window in ms (apply_epoch). If no types specified by GUI, use all current.
if isfield(opt,'apply_epoch') && opt.apply_epoch && isfield(opt,'epoch_window') && numel(opt.epoch_window)>=2 ...
        && isfield(EEG,'event') && ~isempty(EEG.event)
    t_ms = double(opt.epoch_window(:))';
    if numel(t_ms)>=2 && isfinite(t_ms(1)) && isfinite(t_ms(2)) && t_ms(2)>t_ms(1)
        % fprintf('Epoching around %d event type(s), window [%g %g] ms\n', numel(evtTypes), t_ms(1), t_ms(2));
        EEG = pop_epoch(EEG, {}, t_ms/1000, 'epochinfo','yes', 'newname','iEEGLAB epochs');
        EEG = eeg_checkset(EEG);
    else
        warning('Invalid epoch window; skipping epoching.');
    end
    EEG.ieeglab.opt = rmfield(EEG.ieeglab.opt, 'event_filters');
end


% Apply aCAR (Huang et al., 2024) for CCEP data
if isfield(opt,'apply_acar') && opt.apply_acar
    % figure('color','w')
    respData = nan(length(EEG.times), length(EEG.event));
    for iTrial = 1:length(EEG.event)
        % stimElec = extractBefore(EEG.event(iTrial).type, '-');
        respElec = extractAfter(EEG.event(iTrial).type, '-');
        % stimElecIdx = strcmpi({EEG.chanlocs.labels}, stimElec);
        respElecIdx = strcmpi({EEG.chanlocs.labels}, respElec);
        % subplot(2,1,1); hold on

        % plot(EEG.times, squeeze(EEG.data(stimElecIdx,:,iTrial)),'r')
        % plot(EEG.times, squeeze(EEG.data(respElecIdx,:,iTrial)),'k')
        % legend({'stim. elec' 'meas. elec'}); 
        % title("Before aCAR")
        respData(:,iTrial) = squeeze(EEG.data(respElecIdx,:,iTrial));
    end
    figure('color','w'); hold on
    plot(EEG.times, trimmean(respData,20,2), 'LineWidth',2)

    EEG = ieeglab_car(EEG);

    respData = nan(length(EEG.times), length(EEG.event));
    for iTrial = 1:length(EEG.event)
        respElec = extractAfter(EEG.event(iTrial).type, '-');
        respElecIdx = strcmpi({EEG.chanlocs.labels}, respElec);
        respData(:,iTrial) = squeeze(EEG.data(respElecIdx,:,iTrial));
    end

    plot(EEG.times, trimmean(respData,20,2), 'LineWidth',2)
    legend("Before CAR", "After CAR")
    % for iTrial = 1:length(EEG.event)
    %     stimElec = extractBefore(EEG.event(iTrial).type, '-');
    %     respElec = extractAfter(EEG.event(iTrial).type, '-');
    %     stimElecIdx = strcmpi({EEG.chanlocs.labels}, stimElec);
    %     respElecIdx = strcmpi({EEG.chanlocs.labels}, respElec);
    %     subplot(2,1,2); hold on
    %     % plot(EEG.times, squeeze(EEG.data(stimElecIdx,:,iTrial)),'r')
    %     plot(EEG.times, squeeze(EEG.data(respElecIdx,:,iTrial)),'k')
    %     % legend({'stim. elec' 'meas. elec'}); 
    %     title("After aCAR")
    % end
end

% Baseline correction
if isfield(opt,'apply_baseline') && opt.apply_baseline
    EEG = ieeglab_rm_baseline(EEG);
end

