function EEG = ieeglab_load()

EEG = [];

% gui_load1 to files to load and general parameters
[opt, wasCanceled] = ieeglab_gui_load1();
if wasCanceled || isempty(opt)
    return;  % user aborted, exit gracefully
end

% Pull electrode labels from TSV file
if ~isfield(opt, 'elec_tsv') || ~isempty(opt.elec_tsv)
    elecs = readtable(opt.elec_tsv, 'FileType', 'text', 'Delimiter', '\t');
    opt.elec_labels = elecs.name;

    % Check tsv file contains electrode Cartesian (XYZ) coordinates
    idx_x = strcmpi(fieldnames(elecs), 'x');
    if ~any(idx_x) && ~isempty(elecs{:,idx_x})
        warning("No electrodes Cartesian (XYZ) coordinates found in electrodes table (or it was empty). Cannot load electrode coordinates.")
        return
    end
end

% Pull events from TSV file
% if any(contains(lower(opt.analysis_type), {'event-related' 'event' 'evoked' 'ccep'}))
if ~isfield(opt, 'events_tsv') || ~isempty(opt.events_tsv)
    events = readtable(opt.events_tsv, 'FileType', 'text', 'Delimiter', '\t');
    fprintf("%g total events were succesfully loaded from the TSV file.\n", size(events,1))
    opt.events = events; clear events;
    % end
end

% gui_load2 to select channels and events of interest to load (this is just
%  to speed up data loading here to avoid loading unnecessary heavy
%  portions!)
[opt, wasCanceled] = ieeglab_gui_load2(opt);
if wasCanceled || isempty(opt)
    return;  % user aborted, exit gracefully
end


%% Load data (simple)

[~,~,ext] = fileparts(char(opt.dataset_path));

switch lower(ext)
    case '.mefd'

        % --- Metadata
        disp("Loading .mefd metadata...");
        metadata     = ieeglab_load_mefd(opt.dataset_path);
        fs           = double(metadata.time_series_metadata.section_2.sampling_frequency);
        num_samples  = double(metadata.time_series_metadata.section_2.number_of_samples);
        data_len_min = (num_samples / fs) / 60;
        nChan        = numel(opt.chan_list);

        fprintf('Total data length: %.1f min\n', data_len_min);
        fprintf('Sample rate: %.2f Hz\n', fs);
        fprintf('Channels selected: %d\n', nChan);

        % get epoch sample ranges if user wants to load by events
        if isfield(opt,'load_by_events') && opt.load_by_events
            opt = build_event_epochs(opt, fs, num_samples);
            if isempty(opt.epochs_smp)
                warning('No valid epochs matched/within bounds. Loading whole file.');
                opt.epochs_smp = int64([0, num_samples]);  % half-open
            end
        else
            opt.epochs_smp = int64([0, num_samples]);      % half-open
            fprintf('Whole-file load: %.1f min\n', data_len_min);
        end

        % Load signal
        disp("Loading .mefd data...")
        if isempty(gcp('nocreate')), parpool; end
        tic
        [metadata, signal] = ieeglab_load_mefd(opt.dataset_path, [], opt.chan_list, 'samples', opt.epochs_smp);
        toc

        % If epoched, flatten pages to continuous time:
        [~,~,nEp] = size(signal);
        if nEp > 1
            signal = signal(:,:);  % page-wise concat in time
        end

        % Convert back to continuous if epoched during import
        [nChan, nSamp, nEpoch] = size(signal);
        if nEpoch > 1
            signal = signal(:,:);
        end

        % Convert to EEGLAB format
        EEG = eeg_emptyset;
        [~, EEG.filename] = fileparts(opt.dataset_path);
        EEG.filepath = fileparts(opt.dataset_path);
        EEG.data = signal;
        EEG.nbchan = size(EEG.data,1);
        EEG.pnts = size(EEG.data,2);
        EEG.srate = fs;
        EEG = eeg_checkset(EEG);

        % Channel labels from metadata
        for iChan = 1:EEG.nbchan
            EEG.chanlocs(iChan).labels = metadata.time_series_channels(iChan).name;
        end
        EEG = eeg_checkset(EEG);

        % Electrode XYZ coordinates from tsv file
        if ~isempty(opt.elec_tsv)
            EEG = get_elec_coor(EEG, elecs);
        end

        % Load events from ,tsv into EEGLAB
        if ~isempty(opt.events_tsv)
            EEG = load_events(EEG,opt);
        end
        
        % Clear some variables and store opt in EEGLAB structure
        opt = rmfield(opt, 'elec_labels');
        opt = rmfield(opt, 'chan_idx');
        opt = rmfield(opt, 'chan_list');
        EEG.ieeglab.opt = opt;

        % Final check
        EEG = eeg_checkset(EEG);
        % pop_eegplot(EEG,1,1,1);


    case '.nwb'

        % data = nwbRead(opt.dataset_path, 'ignorecache');
        disp("Loading .nwb data...")
        EEG = pop_nwbimport(opt.dataset_path);
        EEG.srate = round(EEG.srate);

        %%%%%%%%%%%%%%%%%% Check events // elecs match selection from GUI!!   %%%%%%%%%%%%%%%%%% 



    case '.vhdr'
        [filepath,filename,ext] = fileparts(opt.dataset_path);
        EEG = pop_loadbv(filepath, sprintf('%s%s', filename, ext), [], [], true);

        data_len_min = EEG.xmax / 60;
        nChan        = numel(opt.chan_list);
        EEG.srate = round(EEG.srate);

        fprintf('Total data length: %.1f min\n', data_len_min);
        fprintf('Sample rate: %.1f Hz\n', EEG.srate);
        fprintf('Channels selected: %d\n', nChan);
        
        % Check sample rate consistency
        fs = round(trimmean( 1000 ./ diff(EEG.times), 20));
        if fs ~= EEG.srate
            error('Discrepancy between sample rate from metadata and sample rate estimated from timestamps.')
        end
        
        % Get epoch sample ranges if user wants to load by events
        if isfield(opt,'load_by_events') && opt.load_by_events
            opt = build_event_epochs(opt, EEG.srate, EEG.pnts);
            if isempty(opt.epochs_smp)
                warning('No valid epochs matched/within bounds. Loading whole file.');
                opt.epochs_smp = int64([0, EEG.pnts]);  % half-open
            end
        else
            opt.epochs_smp = int64([0, EEG.pnts]);      % half-open
            fprintf('Whole-file load: %.1f min\n', data_len_min);
        end
        
        % Load signal
        disp("Loading .vhdr data...")
        EEG = pop_loadbv(filepath, sprintf('%s%s', filename, ext), [opt.epochs_smp(1) opt.epochs_smp(end)], opt.chan_idx, 0);
        
        % Electrode XYZ coordinates from tsv file
        if ~isempty(opt.elec_tsv)
            EEG = get_elec_coor(EEG, elecs);
        end

        % Load events into EEGLAB
        EEG = load_events(EEG, opt);

        % Clear some variables and store opt in EEGLAB structure
        opt = rmfield(opt, 'elec_labels');
        opt = rmfield(opt, 'chan_idx');
        opt = rmfield(opt, 'chan_list');
        EEG.ieeglab.opt = opt;

        % Final check
        EEG = eeg_checkset(EEG);
        % pop_eegplot(EEG,1,1,1);

        % Final check
        EEG = eeg_checkset(EEG);

    otherwise
        error("Sorry, iEEGLAB plugin only supports .nwb and .mefd data formats. Reach out to us with your data type and we'll add it to the list!")
end
