function opt = build_event_epochs(opt, fs, num_samples)
% BUILD_EVENT_EPOCHS  Build [start, stop) sample ranges and event metadata for event-mode loading.
%
% Inputs
%   opt             struct with fields:
%                     .load_by_events   (logical, must be true)
%                     .events           (table with 'onset' [sec])
%                     .event_field      (char/string)   optional, label column
%                     .event_values     (cellstr)       optional, values or {'All'}
%                     .epoch_window_sec [pre post] sec, post > pre
%   fs              sampling rate (Hz)
%   num_samples     total samples in file (int)
%
% Outputs
%   opt.epochs_smp      [nEpoch x 2] int64, 0-based half-open [start, stop)
%   opt.ev_types        {nEpoch x 1} cellstr, per-epoch label ('' if none)
%   opt.ev_lats         [nEpoch x 1] int64, event latency in samples (1-based) AFTER concatenation
%   opt.nSampEpoch      int64, samples per epoch
%   opt.events          with .events possibly filtered (if values provided)

% Defaults (empty → caller decides fallback)
opt.epochs_smp = int64([]);
opt.ev_types   = {};
opt.ev_lats    = int64([]);
opt.nSampEpoch = int64(0);

% Preconditions
if ~isfield(opt,'load_by_events') || ~opt.load_by_events, return; end
if ~isfield(opt,'events') || isempty(opt.events) || ~istable(opt.events), return; end
if ~ismember('onset', opt.events.Properties.VariableNames), return; end
if ~isfield(opt,'epoch_window_sec') || numel(opt.epoch_window_sec)~=2
    error('build_event_epochs:BadWindow','opt.epoch_window_sec must be [pre post] (sec).');
end

% Optional filtering by values (treat {'All'} as no filtering)
if isfield(opt,'event_field') && ~isempty(opt.event_field) ...
        && ismember(opt.event_field, opt.events.Properties.VariableNames) ...
        && isfield(opt,'event_values') && ~isempty(opt.event_values)

    vals  = opt.event_values;
    isAll = iscell(vals) && numel(vals)==1 && ischar(vals{1}) && strcmpi(vals{1},'All');
    if ~isAll
        keep = ismember(string(opt.events.(opt.event_field)), string(vals));
        opt.events = opt.events(keep,:);
    end
end

% Onsets & window
ons = opt.events.onset;
if ~isnumeric(ons), ons = double(string(ons)); end
ons = double(ons(:));
if isempty(ons) || any(~isfinite(ons)), return; end

pre_sec  = double(opt.epoch_window_sec(1));
post_sec = double(opt.epoch_window_sec(2));
if ~(post_sec > pre_sec)
    error('build_event_epochs:BadWindow','epoch_window_sec must satisfy post > pre.');
end

toSamp0    = @(tsec) int64(round(double(tsec) * fs)); % sec -> 0-based samples
opt.nSampEpoch = int64(round((post_sec - pre_sec) * fs));

% Ranges (0-based, half-open)
start0 = toSamp0(ons + pre_sec);
stop0  = start0 + opt.nSampEpoch;

% In-bounds only
inb = (start0 >= 0) & (stop0 <= int64(num_samples));
if ~any(inb), return; end

opt.epochs_smp = [start0(inb) stop0(inb)];

% Labels (optional)
if isfield(opt,'event_field') && ismember(opt.event_field, opt.events.Properties.VariableNames)
    lab = opt.events.(opt.event_field);
    lab = lab(inb,:);
    if iscell(lab), opt.ev_types = lab(:); else, opt.ev_types = cellstr(string(lab(:))); end
else
    opt.ev_types = repmat({''}, nnz(inb), 1);
end

% Event latencies after concatenation (samples, 1-based)
pre_stim_len_sec = -pre_sec;                              % distance from epoch start to event
ev_offset        = int64(round(pre_stim_len_sec * fs));   % samples within each epoch
nEp              = size(opt.epochs_smp,1);
opt.ev_lats          = (int64(0:nEp-1)' .* opt.nSampEpoch) + ev_offset + 1;
end
