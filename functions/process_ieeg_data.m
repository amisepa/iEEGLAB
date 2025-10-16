function [average_ccep,average_ccep_names,tt,srate,crp_out,single_trials] = ...
    process_ieeg_data(mefd_path, events_table, good_channels, channel_areas, baseline_t, t_win_cod, use_CAR, return_single_trials)
% adatapted from ccepPCC_loadAverageSubset.m to solve the readMef3 issues. 
% 
% Read CCEP data directly from a .mefd directory (no readMef3 dependency).
%
% mefd_path: path to 'sub-01_ses-ieeg01_task-ccep_run-01_ieeg.mefd'
% events_table: table with fields: electrical_stimulation_site, onset (sec), status ('good'/other)
% good_channels: numeric indices into channel order returned by mefd_probe
% channel_areas: vector per channel (CRP done where >0)
% baseline_t: [t1 t2] in seconds within epoch
% t_win_cod: [t1 t2] in seconds for CRP window
% use_CAR: 1 to run ccep_CAR64blocks_percent; 0 otherwise
% varargin{1}: return_single_trials (0/1), only if single stim-pair
%
% Cedric Cannard, Aug 2025

% if nargin < 7 || isempty(use_CAR), use_CAR = 1; end
% if ~isempty(varargin) && ~isempty(varargin{1}), return_single_trials = varargin{1}; end
% 
% % try to init FieldTrip if present
% if exist('ft_defaults','file') == 2
%     try
%         ft_defaults;
%     catch ME
%         % Use identifier + formatted message (what the linter wants)
%         warning('process_ieeg_data:FieldTripInit', ...
%                 'ft_defaults failed (%s): %s', ME.identifier, ME.message);
%     end
% end
% % Sanity check that the I/O funcs are visible:
% assert(exist('ft_read_header','file')==2, ...
%     'FieldTrip not on path: ft_read_header not found. addpath(...) and ft_defaults;');
% assert(exist('ft_read_data','file')==2, ...
%     'FieldTrip not on path: ft_read_data not found. addpath(...) and ft_defaults;');
% 
% % Build FieldTrip's Mayo MEF MEX files (macOS/Linux)
% cd('/Users/cedriccannard/Documents/MATLAB/eeglab/plugins/Fieldtrip-lite250523/external/mayo_mef');
% % compile_mexmef;
% setup_mayo_mex
% 
% 
% % Refresh and verify
% rehash toolboxcache;  ft_defaults;
% ft_hastoolbox('mayo_mef', 1);       % should no longer error
% which mefRead -all                  % should list the built MEX

% % ---- probe the .mefd once (names, srate, order, segment map)
% [chan_info, srate] = mefd_probe(mefd_path);
% channel_names = {chan_info.name}';
% nr_channels   = numel(chan_info);


% ---- pre-allocate outputs
average_ccep       = NaN(nr_channels, max(nPairs,1), nTime);
average_ccep_names = cell(max(nPairs,1),1);
crp_out            = [];
single_trials      = [];

% ---- main per-condition loop (loads epochs directly from .tdat)
for kk = 1:nPairs
    fprintf('Loading data for condition %d/%d\n', kk, nPairs);
    these_epochs = find(stim_pair_nr==kk & epochs_include==1);
    if isempty(these_epochs)
        firstIdx = find(stim_pair_nr==kk,1);
        if ~isempty(firstIdx), average_ccep_names{kk} = stim_pair_name{firstIdx}; end
        continue
    end

    average_ccep_names{kk} = stim_pair_name{these_epochs(1)};

    % sample windows for each event
    onset_s   = events_table.onset(these_epochs);
    start_smp = round( (onset_s - epoch_prestim_length_sec) * srate ) + 1;
    stop_smp  = start_smp + nTime - 1;

    % bounds are validated in reader; we'll drop OOB epochs there
    [signaldata, kept] = mefd_read_epochs(mefd_path, chan_info, start_smp, stop_smp);
    if ~all(kept)
        dropped = sum(~kept);
        warning('%d epoch(s) fell outside available data and were dropped.', dropped);
        if ~any(kept), continue; end
    end

    % exclude stimulated channels from CAR pool
    stim_name  = average_ccep_names{kk};
    stimEl1_nm = extractBefore(stim_name,'-');
    stimEl2_nm = extractAfter( stim_name,'-');
    stimEl1_nr = find(ismember(channel_names, stimEl1_nm));
    stimEl2_nr = find(ismember(channel_names, stimEl2_nm));
    good_channels_car = setdiff(good_channels, [stimEl1_nr stimEl2_nr]);

    % optional adjusted CAR
    if use_CAR==1 && ~isempty(good_channels_car)
        perc_channels = 0.2;
        car_timeint   = [0.015 0.500];
        signaldata = ccep_CAR64blocks_percent(signaldata, tt, good_channels_car, perc_channels, car_timeint);
    end

    % baseline subtraction (median)
    samples_base = find(tt>baseline_t(1) & tt<baseline_t(2));
    data_epoch   = ieeg_baselinesubtract(signaldata, samples_base, 'median'); % [chan x trial x time]

    % CRP (limbic-only example as before)
    for ii = 1:numel(channel_areas)
        if channel_areas(ii)>1 && ismember(ii, good_channels_car)
            crp_out(ii,kk).data = squeeze(data_epoch(ii,:,:));
            crp_out(ii,kk).tt   = tt;
            pick = tt>t_win_cod(1) & tt<t_win_cod(2);
            V    = squeeze(data_epoch(ii,:,pick));
            [crp_parms, crp_projs] = CRP_method(V', tt(pick));
            crp_out(ii,kk).crp_parms = crp_parms;
            crp_out(ii,kk).crp_projs = crp_projs;
        else
            crp_out(ii,kk).data = [];
            crp_out(ii,kk).tt   = [];
            crp_out(ii,kk).crp_parms = [];
            crp_out(ii,kk).crp_projs = [];
        end
    end

    % average across trials
    average_ccep(:,kk,:) = squeeze(nanmean(data_epoch, 2));

    % optionally return single trials
    if return_single_trials==1 && nPairs==1
        single_trials = data_epoch; % [chan x trial x time]
    elseif return_single_trials==1 && nPairs>1
        warning('single_trials requested but >1 stim pair; returning empty.');
        single_trials = [];
    end
end

% ---- NaN stimulated electrodes in the averages
for kk = 1:size(average_ccep,2)
    if isempty(average_ccep_names{kk}), continue; end
    el1 = extractBefore(average_ccep_names{kk},'-');
    el2 = extractAfter( average_ccep_names{kk},'-');
    m1 = ismember(channel_names, el1);
    m2 = ismember(channel_names, el2);
    average_ccep(m1,kk,:) = NaN;
    average_ccep(m2,kk,:) = NaN;
end
end

%% Helper function 1


function [chan_info, srate] = mefd_probe(mefd_path, srate_override)
if nargin < 2, srate_override = []; end

all_dirs = dir(mefd_path);
is_timd  = [all_dirs.isdir] & endsWith(lower({all_dirs.name}), '.timd');
d        = all_dirs(is_timd);
if isempty(d)
    error('No *.timd channel folders found in %s', mefd_path);
end

chan_info = struct('name',{},'folder',{},'segments',{},'segd',{},'bytes_per_sample',{}, ...
                   'vconv',{},'units',{},'n_samples',{},'tmet',struct(),'mode',{});
srate = [];

for i = 1:numel(d)
    cdir = fullfile(mefd_path, d(i).name);
    listing = dir(fullfile(cdir, '*'));
    listing = listing(~[listing.isdir]);
    names   = {listing.name}';

    % get channel label from any "<label>-<seg>.<ext>"
    tok = regexp(names, '^([A-Za-z0-9_]+)-(\d+)\.([A-Za-z0-9]+)$', 'tokens', 'once');
    has = ~cellfun('isempty', tok);
    if ~any(has)
        warning('No segment-like files in %s', cdir);
        continue
    end
    base = tok{find(has,1)};
    ch_label = base{1};

    % detect TDAT/SEGD
    tdat_files = names(endsWith(lower(names), '.tdat'));
    segd_files = names(endsWith(lower(names), '.segd'));

    % optional per-channel tmet
    tmet_file = fullfile(cdir, [ch_label '-000000.tmet']);
    tmet = struct(); bps = []; vconv = 1.0; units = '';
    if exist(tmet_file,'file')==2
        txt = fileread(tmet_file);
        tmet.sampling_frequency        = parse_key_num(txt, {'sampling_frequency','sampling rate','Fs'});
        tmet.voltage_conversion_factor = parse_key_num(txt, {'voltage_conversion_factor','voltage conversion'});
        tmet.bytes_per_sample          = parse_key_num(txt, {'bytes_per_sample','bytes per sample'});
        tmet.units                     = parse_key_str(txt, {'units','unit','signal_units'});
        if ~isempty(tmet.sampling_frequency) && isempty(srate), srate = tmet.sampling_frequency; end
        if ~isempty(tmet.bytes_per_sample), bps = tmet.bytes_per_sample; end
        if ~isempty(tmet.voltage_conversion_factor), vconv = tmet.voltage_conversion_factor; end
        if ~isempty(tmet.units), units = tmet.units; end
    end

    chan_info(i).name             = ch_label;
    chan_info(i).folder           = cdir;
    chan_info(i).segments         = tdat_files(:);
    chan_info(i).segd             = segd_files(:);
    chan_info(i).bytes_per_sample = bps;
    chan_info(i).vconv            = vconv;
    chan_info(i).units            = units;
    chan_info(i).tmet             = tmet;
    chan_info(i).n_samples        = NaN;
    chan_info(i).mode             = ternary(~isempty(tdat_files),'tdat','segd');
end

% 1) explicit override wins
if ~isempty(srate_override)
    srate = srate_override;
end

% 2) BIDS sidecar fallback for Fs (your dataset has these)
if isempty(srate)
    json_sidecar = regexprep(mefd_path, '\.mefd$', '.json', 'ignorecase');
    if exist(json_sidecar,'file')==2
        try
            S = jsondecode(fileread(json_sidecar));
            if isfield(S,'SamplingFrequency'), srate = S.SamplingFrequency; end
        catch
            % ignore JSON parse errors here
        end
    end
end

% 3) FieldTrip fallback (requires mayo_mef for SEGDs)
if isempty(srate) && exist('ft_read_header','file')==2
    hdr = ft_read_header(mefd_path);
    srate = hdr.Fs;
end

if isempty(srate)
    error(['Sampling rate not found in .tmet/.json and FieldTrip header unavailable.\n' ...
           'Install mayo\_mef for FieldTrip, or pass srate\_override to process\_ieeg\_data().']);
end
end

function v = parse_key_num(txt, keys)
v = [];
for k = 1:numel(keys)
    m = regexp(txt, [keys{k} '[^\d\-]*([\-]?\d+(\.\d+)?)'], 'tokens', 'ignorecase');
    if ~isempty(m), v = str2double(m{1}{1}); return; end
end
end
function v = parse_key_str(txt, keys)
v = '';
for k = 1:numel(keys)
    m = regexp(txt, [keys{k} '\s*[:=]\s*([^\r\n]+)'], 'tokens', 'ignorecase');
    if ~isempty(m), v = strtrim(m{1}{1}); return; end
end
end
function out = ternary(cond,a,b), if cond, out=a; else, out=b; end, end



%% helper function 2

function [signaldata, kept] = mefd_read_epochs(mefd_path, chan_info, start_smp, stop_smp)
% Returns [nChan x nKeptTrials x nTime]
nChan = numel(chan_info);
nTr   = numel(start_smp);
nTime = mode(stop_smp - start_smp + 1);

% If any channel uses TDAT mode, read with the old raw reader
if any(strcmp({chan_info.mode},'tdat'))
    [signaldata, kept] = mefd_read_epochs_tdat(chan_info, start_smp, stop_smp, nTime);
    return
end

% if FieldTrip is available, prefer it for SEGDs
if exist('ft_read_data','file')==2 && exist('ft_read_header','file')==2
    hdr = ft_read_header(mefd_path);
    total_samp = hdr.nSamples;
    kept = (start_smp >= 1) & (stop_smp <= total_samp);
    idx_keep = find(kept);

    % map our channel order to FieldTrip labels
    [~, chan_idx] = ismember({chan_info.name}', hdr.label);
    if any(chan_idx==0)
        missing = {chan_info{chan_idx==0}.name}; %#ok<CCAT>
        error('Channels missing in FieldTrip header: %s', strjoin(missing, ', '));
    end

    signaldata = NaN(numel(chan_info), numel(idx_keep), mode(stop_smp-start_smp+1), 'double');
    for t = 1:numel(idx_keep)
        tr = idx_keep(t);
        dat = ft_read_data(mefd_path, 'begsample', start_smp(tr), 'endsample', stop_smp(tr), 'chanindx', chan_idx);
        signaldata(:, t, :) = double(dat);
    end
    return
end

% Try a MEF3 mex (if you have one)
if exist('mef3_read_ts_data','file')==3
    kept = true(1,nTr);
    signaldata = NaN(nChan, nTr, nTime, 'double');
    for c = 1:nChan
        for t = 1:nTr
            [x, ok] = mef3_read_ts_data(mefd_path, chan_info(c).name, start_smp(t), stop_smp(t));
            if ~ok, kept(t) = false; continue; end
            signaldata(c,t,:) = double(x); % assume volts; scale if needed
        end
    end
    signaldata = signaldata(:, kept, :);
    return
end

error(['Your .mefd uses SEGDs (no TDAT files), so a MEF3 parser is required.\n' ...
       'Please either:\n' ...
       '  (a) add FieldTrip to your path (ft_read_data), or\n' ...
       '  (b) install/build the MEF3 MATLAB mex (e.g., mef3_read_ts_data).']);
end
