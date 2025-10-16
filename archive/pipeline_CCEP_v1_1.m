%% iEEG cortico-cortical evoked potentials (CCEP) with single pulse stimulation
%  PIPELINE v1.1
%
% Original code: pcc_Allstats_CRP.m (HAPwave)
%
% CODE DEPENDENCIES:
%   - Fieldtrip EEGLAB plugin (to import .mefd data)
%   = Fieldtrip MAYO MEF plugin: https://github.com/jiecui/mef_reader_fieldtrip
%   - HAPwave (https://github.com/MultimodalNeuroimagingLab/HAPwave)
%   - matmef (https://github.com/MaxvandenBoom/matmef)
%
% Cedric Cannard, July 2025

clear; close all; clc

data_path = '/Users/cedriccannard/Downloads/ds004696';
mainDir = '/Users/cedriccannard/Documents/MATLAB/iEEGLAB';

% Add path to iEEGLAB plugin and its subfolders
addpath(genpath(mainDir))
cd(mainDir)

% eeglab; close

% Compile mex MEF (only 1st time) to be able to read MEFD files
make_mex_mef  % original version

%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_subjects    = {'01','02','03','04','05','06','07','08'};    % List of subjects
all_hemi        = {'r','r','r','l','r','l','l','r'};            % List of hemispheres
all_runs        = {'01','01','01','01','01','01','01','01'};    % List of runs
baseline_t      = [-0.5 -0.05];
t_win_crp       = [0 1.5];
stim_currents   = {'4.0 mA', '6.0 mA'};
perform_CAR     = true;                 % perform common average reference (CAR; true) or not (false)
car_chans_ratio = 0.2;                  % portion of channels for CAR  (default = 0.2)    
car_time_int    = [0.015 0.500];       % car time interval (default = [0.015 0.500])
% return_single_trials = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for iSub = 1:length(all_subjects)
iSub = 1;

bids_sub        = all_subjects{iSub};
bids_ses        = '01';
bids_task       = 'ccep';
bids_run        = all_runs{iSub};

% subject folder and file names (events, channels, electrodes)
sub_folder = fullfile(data_path, sprintf('sub-%s',bids_sub), sprintf('ses-ieeg%s', bids_ses),'ieeg');
subject_files = {dir(sub_folder).name}';
tsv_file_channels = fullfile(sub_folder, subject_files{contains(subject_files, 'channels.tsv')});
tsv_file_electrodes = fullfile(sub_folder, subject_files{contains(subject_files, 'electrodes.tsv')});
tsv_file_events = fullfile(sub_folder, subject_files{contains(subject_files, 'events.tsv')});

% Read .TSV data
channels_table   = readtable(tsv_file_channels, 'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
electrodes_table = readtable(tsv_file_electrodes, 'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});

% % Metadata
% [metadata]      = readMef3(fileName); % metadata table
% fs           = metadata.time_series_metadata.section_2.sampling_frequency; % sampling frequency
% nChan     = length(metadata.time_series_channels);                % number of channels
% list of channel names
% channel_names   = cell(nChan,1);
% for kk = 1:nChan
%     channel_names{kk} = metadata.time_series_channels(kk).name;
% end
disp("Extracting metadata from subject json file.")
metadata_file = sprintf('sub-%s_ses-ieeg%s_task-%s_run-%s_ieeg.json', bids_sub, bids_ses, bids_task, bids_run);
info = parse_ieeg_json(fullfile(sub_folder, metadata_file));
fs = info.SamplingFrequency;

% get channel names from BIDS .tsv file
channel_names_tsv  = channels_table.name;

% Pull channels from the .mefd data directory (Cedric bypass)
mefd_folder = fullfile(sub_folder, sprintf('sub-%s_ses-ieeg%s_task-%s_run-%s_ieeg.mefd', bids_sub, bids_ses, bids_task, bids_run));
channel_names_mefd = {dir(mefd_folder).name}';
channel_names_mefd = extractBefore(channel_names_mefd, '.');
channel_names_mefd(ismissing(channel_names_mefd)) = [];
if length(channel_names_mefd) ~= length(channel_names_tsv)
    warning("Different number of channels in .tsv table and in subject data (.mefd) folder")
    % fprintf('Number of channels in .tsc table: %g\n', length(channel_names_tsv))
    % fprintf('Number of channels in .mefd folder: %g\n', length(channel_names_mefd))
    if length(channel_names_mefd) > length(channel_names_tsv)
        fprintf('.mefd folder has %g extra channels relative to .tsv table.\n', length(channel_names_mefd) - length(channel_names_tsv))
    else
        fprintf('.tsv table has %g extra channels relative to .mefd table.\n', length(channel_names_mefd) - length(channel_names_tsv))
    end
else
    fprintf('Same number of sEEG channels detected in .tsv table and .mefd data folder: %g channels.\n', length(channel_names_mefd));
end

% Get channel areas
areas_interest = [17 18 11106 11107 11108 11109 11110 11123 10 53 54 12106 12107 12108 12109 12110 12123 49]; % Destrieux_label
channel_areas  = zeros(size(channel_names_tsv));
for iPair = 1:length(channel_areas)
    % which number has this channel in the electrodes.tsv file?
    thisElPos = find(ismember(electrodes_table.name,channel_names_tsv{iPair}));
    if ismember(electrodes_table.Destrieux_label(thisElPos),areas_interest)
        channel_areas(iPair) = areas_interest(ismember(areas_interest,electrodes_table.Destrieux_label(thisElPos)));
    end
end

% Keep only good sEEG channels --> REMOVING THEM NOW LEADS TO ERROR WITH
% PAIRS LATER OTHERWISE (FIND CORRESPONDING ONES IN PAIRS?)
idx_good = ismember(channels_table.type,{'SEEG'}) & ismember(channels_table.status,'good');
if sum(~idx_good)>0
    bad_channels = channel_names_tsv(~idx_good);
    bad_channels = bad_channels(~cellfun('isempty', bad_channels)); % for when there are empty cells
    warning("%g/%g (%.1f%%) of sEEG channels are labeled as bad.", sum(~idx_good), length(channel_names_tsv), (sum(~idx_good)/length(channel_names_tsv))*100)
    fprintf("Bad channels' labels: %s\n", strjoin(bad_channels, ', '));
    % channel_names_tsv(~idx_good) = [];
    % channel_areas(~idx_good) = [];
    % if sum(idx_good) ~= length(channel_names_tsv)
    %     warning("channel_names and good_channels have different length. Cedric metadata bypass may fail or be unreliable.")
    % end
end

% load event data, keep only events with stim current == 4.0 or 6.0 mA
events_table = readtable(tsv_file_events, 'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
fprintf('Number of events at import: %g\n', size(events_table,1))
fprintf("Keeping only events with stimulation current set to %s.\n", strjoin(stim_currents, ', '))
events_table = bids_clipEvents(events_table,'electrical_stimulation_current', stim_currents);
fprintf('Number of events remaining: %g\n', size(events_table,1))

% Get areas for each stim pair
ccep_stim_names     = unique(events_table.electrical_stimulation_site);
ccep_stim_areas     = zeros(length(ccep_stim_names),2);
for iPair = 1:length(ccep_stim_names)
    if sum(ismember(ccep_stim_names{iPair},'-')) == 1 %
        ccep_stim_areas(iPair,1)   = channel_areas(ismember(channel_names_tsv, extractBefore(ccep_stim_names{iPair},'-')));
        ccep_stim_areas(iPair,2)   = channel_areas(ismember(channel_names_tsv, extractAfter(ccep_stim_names{iPair},'-')));
    else % assume the second - if there is a - in the channel name
        dash_in_name            = find(ismember(ccep_stim_names{iPair},'-'));
        el1_name                = ccep_stim_names{iPair}(1:dash_in_name(2)-1);
        el2_name                = ccep_stim_names{iPair}(dash_in_name(2)+1:end);
        ccep_stim_areas(iPair,1)   = channel_areas(ismember(channel_names_tsv,el1_name));
        ccep_stim_areas(iPair,2)   = channel_areas(ismember(channel_names_tsv,el2_name));
    end
end

% Only keep events while stimulated the limbic network
disp("Keeping only events involving the limbic network.")
ccep_stim_areas_limbic  = ccep_stim_areas(sum(ccep_stim_areas,2)>0,:);
ccep_stim_names_limbic  = ccep_stim_names(sum(ccep_stim_areas,2)>0,:);
events_table = bids_clipEvents(events_table, 'electrical_stimulation_site', ...
    ccep_stim_names_limbic); % include only limbic sites
fprintf('Number of events remaining: %g\n', size(events_table,1))

% Remove events if less than 4 trials of the same type
disp("excluding events with less than 4 trials")
remove_stim_pair        = {};
[uni_stim_sets, bb, cc]   = unique(events_table.electrical_stimulation_site);
for iPair = 1:length(uni_stim_sets)
    if length(find(ismember(events_table.electrical_stimulation_site,uni_stim_sets{iPair})))<4
        disp(['less then 4 events for sub ' int2str(iSub) ' ' uni_stim_sets{iPair}])
        remove_stim_pair =[remove_stim_pair uni_stim_sets(iPair)];
    end
end
if ~isempty(remove_stim_pair)
    warning("Removing %g stim pairs that have less than 4 trials of same type.", length(remove_stim_pair))
    events_table(ismember(events_table.electrical_stimulation_site,remove_stim_pair),:) = [];
end
fprintf('Number of events remaining: %g\n', size(events_table,1))


%% ----------------------- EXTRACT ccepPCC_loadAverageSubset.m ---------%
%  ADAPTEDD FROM ccepPCC_loadAverageSubset.m
% 
% % Only calculate stats for channel_areas>0
% [average_ccep,average_ccep_names,tt,fs,crp_out] = ccepPCC_loadAverageSubset(fileName, ...
%   events_table_clipped,good_channels, channel_areas, baseline_t, ...
%   t_win_crp, use_CAR, return_single_trials);

disp("Initiating preprocessings of sEEG data...")

% Build condition map (canonicalize el1-el2)
stim_pair_nr   = NaN(height(events_table),1);
stim_pair_name = cell(height(events_table),1);
stimEl1 = extractBefore(events_table.electrical_stimulation_site,'-');
stimEl2 = extractAfter( events_table.electrical_stimulation_site,'-');
cond_id = 0;
for k = 1:height(events_table)
    el1 = stimEl1{k}; el2 = stimEl2{k};
    if ~isempty(el1) && ~isempty(el2)
        if ~any(strcmp(stim_pair_name,[el1 '-' el2])) && ~any(strcmp(stim_pair_name,[el2 '-' el1]))
            cond_id = cond_id + 1;
            m = (strcmp(stimEl1,el1) & strcmp(stimEl2,el2)) | (strcmp(stimEl1,el2) & strcmp(stimEl2,el1));
            idx  = find(m);
            for ii = 1:numel(idx), stim_pair_name{idx(ii)} = [el1 '-' el2]; end
            stim_pair_nr(m) = cond_id;
        end
    end
end
nPairs = max([stim_pair_nr; 0]);
fprintf('Number of stim pairs: %g\n', nPairs)

% Epoch parameters (LATENCY MISMATCH WITH OTHER PARAMETERS AT THE TOP???)
epochs_include = ones(size(stim_pair_nr));
epochs_include(~ismember(events_table.status, 'good')) = 0;
epoch_length_sec = 5;     % -2 .. +3 seconds by default
epoch_prestim_length_sec = 2;
nTime = round(epoch_length_sec * fs);
timevec = (0:nTime-1)/fs - epoch_prestim_length_sec;

% Main per-condition loop (loads epochs directly from .tdat)
nChan              = length(channel_names_tsv);
average_ccep       = NaN(nChan, max(nPairs,1), nTime);
average_ccep_names = cell(max(nPairs,1),1);
crp_out            = [];
% single_trials      = [];
% for iPair = 1:nPairs
iPair = 1;

fprintf('Loading data for condition %d/%d\n', iPair, nPairs);

these_epochs = find(stim_pair_nr==iPair & epochs_include==1);
if isempty(these_epochs)
    firstIdx = find(stim_pair_nr==iPair,1);
    if ~isempty(firstIdx)
        average_ccep_names{iPair} = stim_pair_name{firstIdx}; 
    end
    % continue
end

average_ccep_names{iPair} = stim_pair_name{these_epochs(1)};

% sample windows for each event
onset_s   = events_table.onset(these_epochs);
start_smp = round( (onset_s - epoch_prestim_length_sec) * fs ) + 1;
stop_smp  = start_smp + nTime - 1;

%%%%%%%%%%%%%%%%%% IMPORT SEEG SIGNAL FROM MEFD DATA %%%%%%%%%%%%%%%%
% bounds are validated in reader; we'll drop OOB epochs there
% [signal, kept] = mefd_read_epochs(info, start_smp, stop_smp);
disp("Loading .mefd signals for preprocessing...")
[signal, kept] = read_mefd_data(mefd_folder, channel_names_mefd, start_smp, stop_smp);

if ~all(kept)
    dropped = sum(~kept);
    warning('%d epoch(s) fell outside available data and were dropped.', dropped);
    % if ~any(kept), continue; end
end

% Ensure signal has right dimensions
if length(timevec) ~= size(signal, 2) && length(timevec) == size(signal, 3)
    % Assume input is channels × epochs × time → permute
    signal = permute(signal, [1 3 2]);
end

% Check effective data rank at import
dataRank = sum(eig(cov(double(signal(:,:)'))) > 1E-7);  % epoched
% dataRank = sum(eig(cov(double(signal')))>1E-7);    % continuous 
if dataRank < size(signal,1)
    warning("Effective data rank is defficient upon importation : rank = %g and should be %g!", dataRank, size(signal,1))
else
    fprintf("Data are full rank upon importation: %g/%g \n", dataRank, size(signal,1))
end



%%%%%%%%%%%%%%%%%% FILTER %%%%%%%%%%%%%%%%

% highpass filter
cutoff_hp = 0.5;
usefft = true;
signal = filter_ieeg(signal, fs, 'locutoff', cutoff_hp, ...
    'usefftfilt', usefft, 'plotfreqz', false);

% % lowpass filter
% cutoff_lp = fs / 2;
% usefft = false;
% signal_lp = filter_ieeg(signal, fs, 'hicutoff', cutoff_lp, ...
%     'usefftfilt', usefft, 'plotfreqz', false);

% notch filter to remove power line noise
usefft = false;
signal = filter_ieeg(signal, fs, 'locutoff', 59, 'hicutoff', 61,  ...
        'revfilt', 1, 'usefftfilt', usefft, 'plotfreqz', false);


%%%%%%%%%%%%%%%%%% PERFORM CAR %%%%%%%%%%%%%%%%%%
% (excludes stimulated channels from CAR pool)
if perform_CAR

    good_channels = find(idx_good);

    stim_name  = average_ccep_names{iPair};
    stimEl1_nm = extractBefore(stim_name,'-');
    stimEl2_nm = extractAfter( stim_name,'-');
    stimEl1_nr = find(ismember(channel_names_tsv, stimEl1_nm));
    stimEl2_nr = find(ismember(channel_names_tsv, stimEl2_nm));
    good_chans_car = setdiff(good_channels, [stimEl1_nr stimEl2_nr]);

    % ENSURE car_time_int MATCHES EPOCHING ABOVE

    car_preserve_rank = true;      % preserve data rank %%% ALWAYS LOSE 3 RANKS ON 1st TIME BECAUSE OF BLOCK DESIGN DESPITE TRYING RANK PRESERVATION
    [signal, car_info] = apply_ieeg_car(signal, timevec, good_chans_car, car_chans_ratio, car_time_int);
    dataRank = sum(eig(cov(double(signal(:,:)'))) > 1E-7);  % check after CAR
end

% %%%%%%%%%%%%%%%%%% Remove bad trials %%%%%%%%%%%%%%%%%%
% stim_currents = {'6.0 mA', '6.0 mA', '6.0 mA', '6.0 mA', ...
%                  '6.0 mA', '6.0 mA', '6.0 mA', '6.0 mA'}; % -->  get value for each trial from events_table.electrical_stimulation_current
% k_thresh = 6;    % k-threshold multiplier for MAD (default = 6)
% signal = reject_bad_trials_ccep(signal, stim_currents, 1:size(signal,1), k_thresh);

%%%%%%%%%%%%%%%%%% Baseline subtraction %%%%%%%%%%%%%%%%%%
bsl_method = 'median';  % 'median' (default), 'mean', 'trimmed_mean', or '1/f'
bsl_idx = find(timevec>baseline_t(1) & timevec<baseline_t(2));
signal = rm_bsl_ieeg(signal, bsl_idx, bsl_method, fs); % [chan x time]


%%%%%%%%%%%%%%%%%% PLOT CCEP %%%%%%%%%%%%%%%%%%

% Plot all channels as heatmap
plot_ccep(signal, timevec, channel_names_mefd, 'all', [], 0.20);

% Plot single channel
chan_to_plot = 118;
plot_ccep(signal, timevec, channel_names_mefd, 'single', chan_to_plot, []);



%%
%%%%%%%%%%%%%%%%%% CRP %%%%%%%%%%%%%%%%%%
% (limbic-only example as before)
for ii = 1:numel(channel_areas)
    if channel_areas(ii)>1 && ismember(ii, good_chans_car)
        crp_out(ii,iPair).data = squeeze(single_trials(ii,:,:));
        crp_out(ii,iPair).tt   = timevec;
        pick = timevec>t_win_cod(1) & timevec<t_win_cod(2);
        V    = squeeze(single_trials(ii,:,pick));
        [crp_parms, crp_projs] = CRP_method(V', timevec(pick));
        crp_out(ii,iPair).crp_parms = crp_parms;
        crp_out(ii,iPair).crp_projs = crp_projs;
    else
        crp_out(ii,iPair).data = [];
        crp_out(ii,iPair).tt   = [];
        crp_out(ii,iPair).crp_parms = [];
        crp_out(ii,iPair).crp_projs = [];
    end
end

%%%%%%%%%%%%%%%%%% average across trials %%%%%%%%%%%%%%%%%%
average_ccep(:,iPair,:) = squeeze(trimmean(single_trials, 10, 2));

% % optionally return single trials
% if return_single_trials==1 && nPairs==1
%     single_trials = data_epoch; % [chan x trial x time]
% elseif return_single_trials==1 && nPairs>1
%     warning('single_trials requested but >1 stim pair; returning empty.');
%     single_trials = [];
% end
% end

% ---- NaN stimulated electrodes in the averages
for iPair = 1:size(average_ccep,2)
    if isempty(average_ccep_names{iPair}), continue; end
    el1 = extractBefore(average_ccep_names{iPair},'-');
    el2 = extractAfter( average_ccep_names{iPair},'-');
    m1 = ismember(channel_names_tsv, el1);
    m2 = ismember(channel_names_tsv, el2);
    average_ccep(m1,iPair,:) = NaN;
    average_ccep(m2,iPair,:) = NaN;

end

%-------------------- END OF ccepPCC_loadAverageSubset.m -----------------%

% now get stim areas, incorrect order:
[~,x_ind]               = ismember(average_ccep_names,ccep_stim_names_limbic);    % indiced into stim_names from ccep_names
average_ccep_areas      = ccep_stim_areas_limbic(x_ind,:);
% clear ccep_stim_areas_limbic ccep_stim_names_limbic  % these were only used to clip the events table, not for computations

% % Export
% statsFile = fullfile(data_path,'derivatives','stats',['sub-' bids_sub],...
%     ['sub-' bids_sub '_ses-' bids_ses '_task-' bids_task '_run-' bids_run '_crp.mat']);
% save(statsFile,'average_ccep','average_ccep_names','average_ccep_areas','tt','fs','crp_out','channel_names','channel_areas');
% clear average_ccep average_ccep_names tt fs crp_out
disp(['sub ' bids_sub ' stats done'])




