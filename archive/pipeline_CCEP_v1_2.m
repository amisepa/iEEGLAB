%% iEEG cortico-cortical evoked potentials (CCEP) with single pulse stimulation
%  PIPELINE v1.2 (streamlined + alignment-safe)
%
%  Key changes:
%   1) Build a master_table that aligns TSV ↔︎ MEFD (one source of truth).
%   2) Never drop bad channels from names; use masks when loading/processing.
%   3) Canonicalize stim pairs and filter events ONLY after verifying membership & quality.
%   4) All indices used later (CAR, CRP, plots) are in MEFD order via master_table.
%
%  Outputs used later:
%   - fs, mefd_folder
%   - master_table (name, tsv_idx, mefd_idx, type, status, area, good)
%   - channel_areas               (aligned to master_table order)
%   - events_table (filtered to limbic + currents + ≥4 trials, valid/good pairs only)
%   - stim_pair_nr, stim_pair_name (per-event canonical condition ids)
%   - nPairs
%   - timevec, epoch_length_sec, epoch_prestim_length_sec
%   - good_chans_mefd_idx (MEFD-indexed, excluding bad channels)
%   - channel_names (alias for master_table.name; use for plots, reads, etc.)
%
% Notes: 
%   - Channels refer to the individual sEEG recording contacts. 
%   - Pairs are the stimulation electrode pairs. 2 contacts (e.g., A–B) used 
%       together to deliver a single biphasic pulse; each unique pair 
%       defines a "condition" for epoching/averaging/CRP. 
%   - We canonicalize pairs so A–B and B–A are treated as the same stim 
%       site, exclude those two electrodes from CAR, and then compute 
%       responses on all recording channels for each stim pair.
% 
% Cedric Cannard, July 2025

clear; close all; clc

data_path = '/Users/cedriccannard/Downloads/ds004696';
mainDir   = '/Users/cedriccannard/Documents/MATLAB/iEEGLAB';

addpath(genpath(mainDir)); cd(mainDir)
% eeglab; close

% Compile mex MEF (only 1st time) to be able to read MEFD files
make_mex_mef  % original version

%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_subjects    = {'01','02','03','04','05','06','07','08'};    % List of subjects
all_hemi        = {'r','r','r','l','r','l','l','r'};            % List of hemispheres (unused here; kept for completeness)
all_runs        = {'01','01','01','01','01','01','01','01'};    % List of runs

bsl_method     = '1/f';             % or 'median', 'mean', 'trimmed_mean'
baseline_t      = [-0.5 -0.05];     % seconds (must lie within pre-stim)
t_win_crp       = [0.01 1];          % seconds post-stim (used later)
stim_currents   = {'4.0 mA', '6.0 mA'};

perform_CAR     = true;             % perform common average reference (CAR)
car_chans_ratio = 0.2;              % fraction of channels used in CAR
car_time_int    = [0.015 0.500];    % CAR window after stim (sec)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------- Subject selection --------
iSub      = 1;
bids_sub  = all_subjects{iSub};
bids_ses  = '01';
bids_task = 'ccep';
bids_run  = all_runs{iSub};

% -------- BIDS paths / files --------
sub_folder        = fullfile(data_path, sprintf('sub-%s',bids_sub), sprintf('ses-ieeg%s', bids_ses),'ieeg');
subject_files     = {dir(sub_folder).name}';
tsv_file_channels = fullfile(sub_folder, subject_files{contains(subject_files, 'channels.tsv')});
tsv_file_elec     = fullfile(sub_folder, subject_files{contains(subject_files, 'electrodes.tsv')});
tsv_file_events   = fullfile(sub_folder, subject_files{contains(subject_files, 'events.tsv')});
metadata_file     = sprintf('sub-%s_ses-ieeg%s_task-%s_run-%s_ieeg.json', bids_sub, bids_ses, bids_task, bids_run);

% -------- Read metadata & TSVs --------
disp("Reading BIDS TSV/JSON metadata...")
channels_table   = readtable(tsv_file_channels, 'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
electrodes_table = readtable(tsv_file_elec,    'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
events_table     = readtable(tsv_file_events,  'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
info             = parse_ieeg_json(fullfile(sub_folder, metadata_file));
fs               = info.SamplingFrequency;

% -------- MEFD folder & channel list --------
mefd_folder       = fullfile(sub_folder, sprintf('sub-%s_ses-ieeg%s_task-%s_run-%s_ieeg.mefd', bids_sub, bids_ses, bids_task, bids_run));
channel_names_tsv = channels_table.name;
channel_names_mef = {dir(mefd_folder).name}';
channel_names_mef = extractBefore(channel_names_mef, '.');                % strip extensions
channel_names_mef(ismissing(channel_names_mef)) = [];

% -------- Align TSV ↔︎ MEFD into a master table --------
% Stable intersection keeps TSV order; unmatched channels are dropped with warning.
[common_chans, idx_tsv, idx_mef] = intersect(channel_names_tsv, channel_names_mef, 'stable');
if numel(common_chans) ~= numel(channel_names_tsv)
    warn_drop = setdiff(channel_names_tsv, common_chans, 'stable');
    if ~isempty(warn_drop)
        warning("Dropping %d TSV channels not found in MEFD: %s", numel(warn_drop), strjoin(warn_drop(~cellfun('isempty',warn_drop))', ', '));
    end
end
extra_mef = setdiff(channel_names_mef, common_chans, 'stable');
if ~isempty(extra_mef)
    warning("Ignoring %d MEFD channels absent from TSV: %s", numel(extra_mef), strjoin(extra_mef(~cellfun('isempty',extra_mef))', ', '));
end

master_table = table;
master_table.name     = channel_names_tsv(idx_tsv);
master_table.tsv_idx  = idx_tsv;
master_table.mefd_idx = idx_mef;

% Attach channel metadata aligned to master_table
master_table.type   = channels_table.type(idx_tsv);
% Some datasets may lack 'status'; treat missing as 'good'
if ismember('status', channels_table.Properties.VariableNames)
    master_table.status = channels_table.status(idx_tsv);
else
    master_table.status = repmat("good", height(master_table), 1);
end

% Map anatomical area (Destrieux_label) to master_table rows
% (ensure electrodes_table is filtered to rows present in master_table)
elec_map = containers.Map(electrodes_table.name, electrodes_table.Destrieux_label);
area_vec = zeros(height(master_table),1);
for iChan = 1:height(master_table)
    nm = master_table.name{iChan};
    if isKey(elec_map, nm)
        v = elec_map(nm);
        if ~ismissing(v) && ~isnan(v), area_vec(iChan) = v; end
    end
end
master_table.area = area_vec;

% Define good channel mask (keep all names; mask used at processing time)
master_table.good = (master_table.type == "SEEG") & (master_table.status == "good");
good_chans_mefd_idx = master_table.mefd_idx(master_table.good);
channel_names       = master_table.name;          % <- use this everywhere for plotting/reads
channel_areas       = master_table.area;          % <- aligned to master_table rows

% -------- Event filtering: currents → limbic → valid & good pairs → n≥4 --------
fprintf('Events at import: %d\n', height(events_table))

% 1) Keep only specific stimulation currents
fprintf("Filtering events by current: %s\n", strjoin(stim_currents, ', '))
if ismember('electrical_stimulation_current', events_table.Properties.VariableNames)
    events_table = bids_clipEvents(events_table, 'electrical_stimulation_current', stim_currents);
else
    warning("Column 'electrical_stimulation_current' not found; skipping current filter.");
end
fprintf('Events after current filter: %d\n', height(events_table))

% 2) Canonicalize pair strings and verify channel presence/quality
%    - Canonical form: "A-B" with A,B as they appear; A-B and B-A treated as the same condition.
stim_site = events_table.electrical_stimulation_site;
stimEl1   = extractBefore(stim_site, '-');
stimEl2   = extractAfter( stim_site, '-');

is_valid_pair = false(height(events_table),1);
is_good_pair  = false(height(events_table),1);
canon_pair    = strings(height(events_table),1);

name_to_row = containers.Map(master_table.name, num2cell(1:height(master_table)));  % master index

for k = 1:height(events_table)
    a = strtrim(stimEl1{k});
    b = strtrim(stimEl2{k});

    ok = ~isempty(a) && ~isempty(b) && isKey(name_to_row, a) && isKey(name_to_row, b);
    is_valid_pair(k) = ok;
    if ~ok
        canon_pair(k) = "";
        continue
    end

    ia = name_to_row(a);
    ib = name_to_row(b);
    is_good_pair(k) = master_table.good(ia) && master_table.good(ib);

    % Canonicalize (unordered): lower index first → stable across the pipeline
    if ia <= ib
        canon_pair(k) = string(master_table.name{ia}) + "-" + string(master_table.name{ib});
    else
        canon_pair(k) = string(master_table.name{ib}) + "-" + string(master_table.name{ia});
    end
end
events_table.canon_pair   = canon_pair;
events_table.is_valid_pair = is_valid_pair;
events_table.is_good_pair  = is_good_pair;

% 3) Limbic-only: keep events where both electrodes have area ∈ areas_interest
areas_interest = [17 18 11106 11107 11108 11109 11110 11123 10 53 54 12106 12107 12108 12109 12110 12123 49]; % Destrieux_label
in_limbic = false(height(events_table),1);
for k = 1:height(events_table)
    if ~events_table.is_valid_pair(k), continue; end
    a = stimEl1{k}; b = stimEl2{k};
    ia = name_to_row(a); ib = name_to_row(b);
    in_limbic(k) = ismember(master_table.area(ia), areas_interest) || ismember(master_table.area(ib), areas_interest);
end

events_table = events_table(events_table.is_valid_pair & events_table.is_good_pair & in_limbic, :);
fprintf('Events after limbic + valid/good pair filter: %d\n', height(events_table))


% 4) Drop stim pairs with < 4 trials
if height(events_table)>0
    [uniq_pairs,~,pair_ids] = unique(events_table.canon_pair);
    grp_counts = accumarray(pair_ids, 1);
    keep_pairs = uniq_pairs(grp_counts >= 4);
    events_table = events_table(ismember(events_table.canon_pair, keep_pairs), :);
end
fprintf('Events after ≥4-trials filter: %d\n', height(events_table))

% -------- Build condition map (per-event) --------
% Assign integer condition IDs for each canonical pair (unordered pairs merged).
stim_pair_name = cell(height(events_table),1);
stim_pair_nr   = NaN(height(events_table),1);
cond_id = 0;
for k = 1:height(events_table)
    cp = events_table.canon_pair{k};
    if cp == "", continue; end
    if ~any(strcmp(stim_pair_name, cp))
        cond_id = cond_id + 1;
        stim_pair_name(strcmp(events_table.canon_pair, cp)) = {cp};
        stim_pair_nr(strcmp(events_table.canon_pair, cp))   = cond_id;
    end
end
nPairs = max([stim_pair_nr; 0]);
fprintf('Number of stim pairs: %d\n', nPairs)

% -------- Epoching parameters --------
% Choose epoch windows up-front and verify baseline lies inside pre-stim.
epoch_length_sec         = 5.0;    % total window length (sec)
epoch_prestim_length_sec = 2.0;    % seconds before stim (timevec starts at -2)
nTime   = round(epoch_length_sec * fs);
timevec = (0:nTime-1)/fs - epoch_prestim_length_sec;

if ~(baseline_t(1) >= -epoch_prestim_length_sec && baseline_t(2) <= 0)
    warning('baseline_t [%g %g] must lie within [-%g, 0].', baseline_t(1), baseline_t(2), epoch_prestim_length_sec);
end
if ~(car_time_int(1) >= 0 && car_time_int(2) <= (epoch_length_sec - epoch_prestim_length_sec))
    warning('car_time_int [%g %g] should lie within [0, %g].', car_time_int(1), car_time_int(2), epoch_length_sec-epoch_prestim_length_sec);
end

% Mark includable epochs (require events "status==good" if present)
if ismember('status', events_table.Properties.VariableNames)
    st = string(events_table.status);     % robust for char/categorical
    epochs_include = (st == "good");      % logical
    epochs_include(ismissing(epochs_include)) = false;   % treat missing as not-good
else
    epochs_include = true(height(events_table),1);       % if no status column, include all
end


% -------- Handy lookups for later stages --------
% MEFD indices of stim channels for each condition (for CAR exclusions, NaN in averages, etc.)
% Build per-pair mappings in master-table (row) space, plus canonical names.
pair_to_rows  = cell(nPairs,1);   % [rowA, rowB] in master_table/channel order
pair_to_mefd  = cell(nPairs,1);   % [mefd_idx_A, mefd_idx_B] (for record only)
pair_to_names = cell(nPairs,1);   % {'A-B'} canonical pair name
for p = 1:nPairs
    rows = find(stim_pair_nr == p);
    if isempty(rows), continue; end

    % Use the canonical name we already wrote into events_table.canon_pair
    cp  = events_table.canon_pair{rows(1)};
    if cp == "", continue; end

    el1 = extractBefore(cp,'-');
    el2 = extractAfter( cp,'-');

    % Master-table ROW indices (these index directly into signal/average_ccep)
    ia  = name_to_row(el1);
    ib  = name_to_row(el2);
    pair_to_rows{p} = [ia, ib];

    % MEFD indices (kept just for provenance; do NOT use as array indices)
    pair_to_mefd{p} = [master_table.mefd_idx(ia), master_table.mefd_idx(ib)];

    % Canonical names (optional but handy downstream)
    pair_to_names{p} = cp;
end

% If you want average_ccep_names to be fully populated up-front:
average_ccep_names = pair_to_names(:);

disp("Setup complete: master_table aligned, events filtered, conditions built, epoching defined.")



%% ----------------------- CCEP extraction & preprocessing -----------------------
% Requires from previous section:
% - fs, mefd_folder, master_table (fields: name, good), name_to_row (map)
% - events_table (with canon_pair), stim_pair_nr, nPairs
% - epoch_length_sec, epoch_prestim_length_sec, timevec
% - perform_CAR, car_chans_ratio, car_time_int
% - baseline_t, t_win_crp

disp("Initiating preprocessing of sEEG data...")

channel_names = master_table.name;          % MEFD-aligned order
nChan         = height(master_table);
nTime         = numel(timevec);

% Preallocate outputs
average_ccep        = NaN(nChan, max(nPairs,1), nTime);
average_ccep_names  = cell(max(nPairs,1),1);

% Helper: baseline index once
bsl_idx = find(timevec > baseline_t(1) & timevec < baseline_t(2));
if isempty(bsl_idx)
    warning('Baseline window [%g %g]s produced empty index. Check epoch window.', baseline_t(1), baseline_t(2));
end

% Loop conditions
% for iPair = 1:max(nPairs,1)
iPair = 8;

fprintf('Loading data for condition %d/%d\n', iPair, max(nPairs,1));

% ---- select epochs for this condition ----
these_epochs = find((stim_pair_nr == iPair));
if isempty(these_epochs)
    % If no epochs survived earlier filtering, skip
    firstIdx = find(stim_pair_nr == iPair, 1);
    if ~isempty(firstIdx)
        average_ccep_names{iPair} = events_table.canon_pair{firstIdx};
    end
    % continue
end

% Respect per-event include mask if present (epochs_include should be logical)
if exist('epochs_include','var') && ~isempty(epochs_include)
    these_epochs = these_epochs(epochs_include(these_epochs));
    % if isempty(these_epochs), continue; end
end

% Canonical name for this pair
average_ccep_names{iPair} = events_table.canon_pair{these_epochs(1)};

% Stim electrodes (row indices in master_table/channel_names order)
el1 = extractBefore(average_ccep_names{iPair},'-');
el2 = extractAfter( average_ccep_names{iPair},'-');
ia  = name_to_row(el1);
ib  = name_to_row(el2);
stim_rows = [ia, ib];

% ---- per-epoch sample windows ----
onset_s   = events_table.onset(these_epochs);
start_smp = round( (onset_s - epoch_prestim_length_sec) * fs ) + 1;
stop_smp  = start_smp + nTime - 1;

% ---- read MEFD epochs (channels follow channel_names order) ----
[signal, kept] = read_mefd_data(mefd_folder, channel_names, start_smp, stop_smp);
if ~all(kept)
    warning('%d/%d epoch(s) fell outside available data and were dropped (pair %d).', sum(~kept), numel(kept), iPair);
end
% Standardize dims to [chan x time x trials]
if ndims(signal) == 2
    % some readers may collapse trial dim if single trial; expand
    signal = reshape(signal, size(signal,1), size(signal,2), 1);
end
% If reader returned [chan x trials x time], permute to [chan x time x trials]
if size(signal,2) ~= nTime && size(signal,3) == nTime
    signal = permute(signal, [1 3 2]);
end
if size(signal,2) ~= nTime
    error('Epoch length mismatch: expected nTime=%d, got %d.', nTime, size(signal,2));
end

% ---- filtering (HP + notch), vectorized over trials ----
% (1) High-pass (FFT for stability on long windows)
cutoff_hp = 0.5;
sig2D = reshape(signal, nChan, []);  % [chan x (time*trials)]
sig2D = filter_ieeg(sig2D, fs, 'locutoff', cutoff_hp, 'usefftfilt', true,  'plotfreqz', false);

% (2) Notch 60 Hz (IIR default)
sig2D = filter_ieeg(sig2D, fs, 'locutoff', 59, 'hicutoff', 61, 'revfilt', 1, 'usefftfilt', false, 'plotfreqz', false);
signal   = reshape(sig2D, nChan, nTime, []);

% ---- CAR (exclude stim channels + bad channels) ----
if perform_CAR
    % Analysis-eligible channels for this pair: good & non-stim
    good_rows      = find(master_table.good);
    stim_rows      = pair_to_rows{iPair};        % <-- row indices, not mefd
    good_chans_car = setdiff(good_rows, stim_rows);


    if isempty(good_chans_car)
        warning('No channels available for CAR pool in pair %d; skipping CAR.', iPair);
    else
        [signal, car_info] = apply_ieeg_car(signal, timevec, good_chans_car, car_chans_ratio, car_time_int);
    end
end

% ---- remove baseline across all trials at once (expects [chan x time x trials]) ----
if ndims(signal) == 2
    % just in case the reader collapsed trial dim for a single kept epoch
    signal = reshape(signal, size(signal,1), size(signal,2), 1);
end
signal = rm_bsl_ieeg(signal, bsl_idx, bsl_method, fs);

% Plot 1: heatmap of the average across trials
plot_ccep(signal, timevec, channel_names, 'all', [], 0.20);

% Plot 2: single-channel overlay of all trials + mean
chan_to_plot = min(118, nChan);   % select channl to plot here
plot_ccep(signal, timevec, channel_names, 'single', chan_to_plot, []);


%% %%%%%%%%%%%%%%%%% Canonical response parameterization (CRP) %%%%%%%%%%%%%%%%%%
% see file kjm_CRP_code_and_data_description.pdf
% (limbic-only example as before)
% Requires: signal [chan x time x trials], timevec, master_table, pair_to_mefd{iPair}, channel_areas

% CRP analysis window (avoid immediate stim artifact)
pick = (timevec >= t_win_crp(1)) & (timevec <= t_win_crp(2));
if ~any(pick)
    warning('CRP window [%g %g] selects no samples; check epoch/timevec.', t_win_crp(1), t_win_crp(2));
end

% Analysis-eligible channels for this pair: good & non-stim
good_rows      = find(master_table.good);
stim_rows      = pair_to_mefd{iPair};        % MEFD indices of stim electrodes
good_chans_car = setdiff(good_rows, stim_rows);

% Initialize container on first pair
if iPair == 1
    crp_out = repmat(struct('data',[],'tt',[],'crp_parms',[],'crp_projs',[]), nChan, nPairs);
end

progressbar('Performing Canonical response parameterization (CRP) on all channels')
for iChan = 1:nChan
    if channel_areas(iChan) ~= 0 && ismember(iChan, good_chans_car)

        % V_full: [time x trials] for this channel
        V_full = squeeze(signal(iChan, :, :));        % [nTime x nTrials]
        if isvector(V_full), V_full = V_full(:); end

        % Store trials x time for inspection (keeps your original habit)
        crp_out(iChan,iPair).data = V_full.';      % [trials x time]
        crp_out(iChan,iPair).tt   = timevec;

        % Guards: enough time samples & trials for CRP
        has_enough_time   = sum(pick) >= 10;
        has_enough_trials = size(V_full,2) >= 2;

        if any(pick) && has_enough_time && has_enough_trials
            V_win = V_full(pick, :);            % [Tsel x K]  (time x trials)
            % [crp_parms, crp_projs] = run_CRP(V_win, timevec(pick));
            [crp_parms, crp_projs] = run_CRP(V_win, timevec(pick), ...
                    struct('label', sprintf('Chan %d, Pair %d', iChan, iPair), 'verbose', true));
            crp_out(iChan,iPair).crp_parms = crp_parms;
            crp_out(iChan,iPair).crp_projs = crp_projs;
        else
            crp_out(iChan,iPair).crp_parms = [];
            crp_out(iChan,iPair).crp_projs = [];
        end
    else
        crp_out(iChan,iPair).data = [];
        crp_out(iChan,iPair).tt   = [];
        crp_out(iChan,iPair).crp_parms = [];
        crp_out(iChan,iPair).crp_projs = [];
    end

    progressbar(iChan / nChan)
end


% Average across trials (kept here for clarity; same as earlier avg_pair)
average_ccep(:, iPair, :) = squeeze(trimmean(signal, 10, 3));   % [chan x time]

% Blank out the two stimulated electrodes for each stim pair in our averaged
% CCEP cube. For each condition, pair_to_rows{iPair} gives the row indices
% of the two stim contacts in channel_names/master_table order.
% We set average_ccep(stim_rows, iPair, :) = NaN so those rows don't show
% huge stimulation artifacts or saturations in heatmaps/stats.
% (This does not affect CRP results; it only masks the averages.)
for iPair = 1:size(average_ccep,2)
    if iPair <= numel(pair_to_rows) && ~isempty(pair_to_rows{iPair})
        stim_rows = pair_to_rows{iPair};     % row indices
        average_ccep(stim_rows, iPair, :) = NaN;
    end
end

% For each stim pair, record the anatomical areas of the two electrodes.
% pair_to_rows{p} gives the row indices of the stim contacts,
% and master_table.area holds their Destrieux labels.
% Result: average_ccep_areas is [nPairs x 2], with one area per contact.
average_ccep_areas = NaN(nPairs, 2);
for p = 1:nPairs
    if p > numel(pair_to_rows) || isempty(pair_to_rows{p}), continue; end
    r = pair_to_rows{p};   % [rowA, rowB]
    average_ccep_areas(p, :) = [master_table.area(r(1)), master_table.area(r(2))];
end


% % Export
% statsFile = fullfile(data_path,'derivatives','stats',['sub-' bids_sub],...
%     ['sub-' bids_sub '_ses-' bids_ses '_task-' bids_task '_run-' bids_run '_crp.mat']);
% save(statsFile,'average_ccep','average_ccep_names','average_ccep_areas','tt','fs','crp_out','channel_names','channel_areas');
% clear average_ccep average_ccep_names tt fs crp_out

% end
