function EEG = load_events(EEG, opt)

EEG.event = [];

% % Compute ev_types / ev_lats
% opt = build_event_epochs(opt, EEG.srate, EEG.pnts);

% % Get epoch sample ranges
% opt = build_event_epochs(opt, fs, num_samples);

% % integrate into EEGLAB dataset
% nEv = numel(opt.ev_lats);
% fprintf("Loading %g events into EEGLAB dataset...\n", nEv);
% for iEv = 1:nEv
%     EEG.event(iEv).type    = opt.ev_types{iEv};
%     EEG.event(iEv).latency = double(opt.ev_lats(iEv));  % samples, 1-based
% end


% whole-file: take from original events
% else
    if isfield(opt,'events') && ~isempty(opt.events) && ismember('onset', opt.events.Properties.VariableNames)
        nEv = height(opt.events);
        fprintf("Loading %g events into EEGLAB dataset...\n", nEv)
        for iEv = 1:nEv
            % Prefer selected field if present; fallback to electrical_stimulation_site
            if isfield(opt,'event_field') && ~isempty(opt.event_field) && ismember(opt.event_field, opt.events.Properties.VariableNames)
                ev_label = opt.events.(opt.event_field){iEv};
            % elseif ismember('electrical_stimulation_site', opt.events.Properties.VariableNames)
            %     val = opt.events.electrical_stimulation_site;
            else
                error("Event field not specified. Please select the event field to use for analysis in the GUI and try again.");
            end

            % Type as string/cellstr
            if iscell(ev_label)
                EEG.event(iEv).type = ev_label;
            else         
                EEG.event(iEv).type = char(string(val(iEv)));
            end

            % Latency in samples (1-based)
            EEG.event(iEv).latency = round(double(opt.events.onset(iEv)) * EEG.srate) + 1;
        end
    end
% end

EEG = eeg_checkset(EEG, 'eventconsistency');