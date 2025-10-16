function EEG = load_events(EEG, opt)

EEG.event = [];

% when epoched+concatenated: use ev_types / ev_lats computed above
if isfield(opt,'load_by_events') && opt.load_by_events
    nEv = numel(opt.ev_lats);
    fprintf("Loading %g events into EEGLAB dataset...\n", nEv);
    for iEv = 1:nEv
        EEG.event(iEv).type    = opt.ev_types{iEv};
        EEG.event(iEv).latency = double(opt.ev_lats(iEv));  % samples, 1-based
    end

% whole-file: take from original events
else
    if isfield(opt,'events') && ~isempty(opt.events) && ismember('onset', opt.events.Properties.VariableNames)
        nEv = height(opt.events);
        fprintf("Loading %g events into EEGLAB dataset...\n", nEv)
        for iEv = 1:nEv
            % Prefer selected field if present; fallback to electrical_stimulation_site
            if isfield(opt,'event_field') && ~isempty(opt.event_field) && ismember(opt.event_field, opt.events.Properties.VariableNames)
                val = opt.events.(opt.event_field);
            elseif ismember('electrical_stimulation_site', opt.events.Properties.VariableNames)
                val = opt.events.electrical_stimulation_site;
            else
                error("Event field not specified and default field 'electrical_stimulation_site' not found. Please select the field in the GUI and try again.");
            end
            % Type as string/cellstr
            if iscell(val)
                EEG.event(iEv).type = val{iEv};
            else         
                EEG.event(iEv).type = char(string(val(iEv)));
            end

            % Latency in samples (1-based)
            EEG.event(iEv).latency = round(double(opt.events.onset(iEv)) * EEG.srate) + 1;
        end
    end
end

EEG = eeg_checkset(EEG, 'eventconsistency');
