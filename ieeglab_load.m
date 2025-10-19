function EEG = ieeglab_load(EEG)

% check if dataset is already epoched
if ndims(EEG.data) == 3
    % epoched = true;
    % if ~isfield(EEG, 'trials')
    %     EEG.trials = size(EEG.data,3);
    % end
else
    % epoched = false;
    if ~isfield(EEG, 'trials') || isempty(EEG.trials) || EEG.trials == 0
        EEG.trials = 1;
    end
end

% gui_load1 to files to load and general parameters
[opt, wasCanceled] = ieeglab_gui_load1(EEG.filepath);
if wasCanceled || isempty(opt)
    return;  % user aborted, exit gracefully
end

% Pull electrode labels from TSV file
if ~isempty(opt.elec_tsv)
    elecs = readtable(opt.elec_tsv, 'FileType', 'text', 'Delimiter', '\t');
    % opt.elecs = elecs;
    opt.elec_labels = elecs.name;
    opt.elec_labels = strtrim(regexprep(opt.elec_labels, '''', '')); % remove apostrophes inside the labels in some datasets

    % Check tsv file contains electrode Cartesian (XYZ) coordinates
    idx_x = strcmpi(fieldnames(elecs), 'x');
    if ~any(idx_x) && ~isempty(elecs{:,idx_x})
        warning("No electrodes Cartesian (XYZ) coordinates found in electrodes table (or it was empty). Cannot load electrode coordinates.")
        return
    end
end

% Pull events from TSV file or EEGLAB dataset
if ~isempty(opt.events_tsv)
    events = readtable(opt.events_tsv, 'FileType', 'text', 'Delimiter', '\t');
    fprintf("%g total events were imported from the .tsv event file.\n", size(events,1))
    opt.events = events; 
else
    disp("You did not select a .tsv event file. Checking if you EEGLAB dataset already contains events...")
    if ~isempty(EEG.event)
        disp("EEGLAB dataset does contain events.")
        opt.events = struct2table(EEG.event);
    else
        warning("No events detected in the EEGLAB dataset you loaded --> Analysis mode: CONTINUOUS.")
    end
end

% Electrode XYZ coordinates from TSV elec data
if ~isempty(opt.elec_tsv)
    EEG = get_elec_coor(EEG, elecs);
    opt.elec_labels = {EEG.chanlocs.labels};
end

% Abort if no channels left
if EEG.nbchan == 0
    error("No electrodes left in dataset.")
end

% Abort if there are no electrode labels anywhere
if ~isfield(EEG.chanlocs, 'X') || isempty([EEG.chanlocs.X])
    error("No electrode locations. Cannot proceed with analysis. You must load a valid .tsv file containing the electrodes' XYZ coordinates.")
end


% gui_load2 to select channels of interest
[opt, wasCanceled] = ieeglab_gui_load2(opt);
if wasCanceled || isempty(opt)
    return  % user aborted, exit gracefully
end


% pull event types and latencies from TSV event data
if ~isempty(opt.events_tsv)
    disp("Loading TSV events into EEGLAB dataset...")
    
    col = opt.event_field;
    ev_types = opt.events.(col);
    ev_lats = opt.events.onset;
    for iEv = 1:length(ev_types)
        if iscell(ev_types(iEv))
            EEG.event(iEv).type = ev_types{iEv};
        else
            EEG.event(iEv).type = ev_types(iEv);
        end
        EEG.event(iEv).latency = ev_lats(iEv) * EEG.srate;
    end
end


% Keep only events of interest (if selected by user)
if isfield(opt, 'event_values') && ~isempty(opt.event_values)
    trials_to_rem = ~contains({EEG.event.type}, opt.event_values);
    % EEG = pop_select(EEG, 'rmtrial', trials_to_rem); % only for epoched data
    EEG.event(trials_to_rem) = [];
    EEG = eeg_checkset(EEG, 'eventconsistency');
    EEG = eeg_checkset(EEG);
end

% Keep only electrodes of interest
if isfield(opt, 'chan_idx') && ~isempty(opt.chan_idx) || ...
       isfield(opt, 'chan_list') && ~isempty(opt.chan_list) && ...
       length(opt.chan_list) ~= EEG.nbchan

    opt.chan_list = cellstr(opt.chan_list);  % ensure cellstr
    % opt.chan_list = opt.chan_list(:)';       % row cell array


    EEG = pop_select(EEG, 'channel', opt.chan_list);
    % EEG = pop_select(EEG, 'channel', opt.chan_list, 'sorttrial', 'off');
    % EEG = pop_select( EEG, 'channel',{'ROC1','ROC2','ROC3','ROC4','ROC5','ROC6','ROC7','ROC8','ROC9','ROC10','ROC11','ROC12','ROC13','ROC14','ROC15'});

    % opt = rmfield(opt, 'chan_idx');
    % opt = rmfield(opt, 'chan_list');
end

% % Clear vars
% opt = rmfield(opt, 'elec_tsv');
% opt = rmfield(opt, 'elecs');
% opt = rmfield(opt, 'elec_labels');
% opt = rmfield(opt, 'events');
% opt = rmfield(opt, 'event_field');
% opt = rmfield(opt, 'event_values');
% opt = rmfield(opt, 'events_tsv');

EEG.ieeglab.opt = opt;

% Final check
EEG = eeg_checkset(EEG);
% pop_eegplot(EEG,1,1,1);
