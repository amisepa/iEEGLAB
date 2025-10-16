function [EEG, out] = ieeglab_car(EEG)
% ieeglab_car  aCAR with per-epoch stim-contact exclusion (C x T x N).
% Uses GUI options in EEG.ieeglab.opt:
%   apply_acar    : logical
%   acar_fraction : default 0.20
%   acar_timewin  : [start end] in ms (default [15 500])
%
% Excludes stim contact(s) per epoch using, in order of precedence:
%   EEG.ieeglab.opt.events.(electrodes_involved_onset | electrical_stimulation_site | stim_electrodes)
%   or EEG.ieeglab.opt.events.(stim_type + stim_num)
%   or equivalent fields on EEG.epoch
%
% Output:
%   EEG.data re-referenced; EEG.ref = 'aCAR'
%   out      per-block diagnostics

% ----- options -----
opt = struct('apply_acar', false, 'acar_fraction', 0.20, 'acar_timewin', [15 500]);
if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'opt') && isstruct(EEG.ieeglab.opt)
    f = EEG.ieeglab.opt;
    if isfield(f,'apply_acar'),    opt.apply_acar = logical(f.apply_acar); end
    if isfield(f,'acar_fraction') && ~isempty(f.acar_fraction)
        opt.acar_fraction = max(0, min(1, double(f.acar_fraction)));
    end
    if isfield(f,'acar_timewin') && numel(f.acar_timewin)>=2
        opt.acar_timewin = double(f.acar_timewin(1:2));
    end
end
if ~opt.apply_acar, out = struct(); return; end

% ----- data shape: C x T x N -----
X = EEG.data;
if ndims(X)==2, X = reshape(X, size(X,1), size(X,2), 1); end
[C,T,N] = size(X);

% ----- time window mask in ms -----
if isfield(EEG,'times') && ~isempty(EEG.times) && numel(EEG.times)==T
    t_ms  = double(EEG.times(:));
    tmask = (t_ms >= opt.acar_timewin(1)) & (t_ms <= opt.acar_timewin(2));
    if ~any(tmask)
        warning('aCAR: No samples in [%g %g] ms; using all samples.', opt.acar_timewin(1), opt.acar_timewin(2));
        tmask = true(T,1);
    end
else
    tmask = true(T,1);
end

% ----- map per-epoch stim exclusions (indices into channels) -----
labels = string({EEG.chanlocs.labels});
stim_excl = local_collect_stim_exclusions(EEG, labels, N); % cell{N} of double indices
stim_union = unique(cat(2, stim_excl{:}));                  % union across epochs

% ----- diagnostics & rank proxy -----
try
    mean_corr = mean(triu(corr(mean(X,3)') ,1),'all','omitnan');
    if mean_corr < 0.01
        warning('aCAR: Very low inter-channel correlation; data may already be referenced.');
    end
catch, end
r_before = local_rank_proxy(X);

% ----- build 64-ch blocks -----
blkSize = 64;
nBlk = ceil(C/blkSize);
out = repmat(struct('channels_set',[],'car_channels',[],'var_sel_thr',NaN,'n_sel',0), nBlk, 1);
for b = 1:nBlk
    idx = (b-1)*blkSize + (1:blkSize);
    out(b).channels_set = idx(idx<=C);
end

% ----- variance-based selector per block (exclude stim_union from pool) ---
frac = opt.acar_fraction;
for b = 1:nBlk
    idxSet = out(b).channels_set;
    if isempty(idxSet), continue; end

    % EXCLUDE stim channels from the selection pool (validated behavior)
    idxPool = setdiff(idxSet, stim_union);
    if isempty(idxPool), idxPool = idxSet; end % fallback

    Xi = X(idxPool, tmask, :);                % [cp x Tw x N]
    Xi = reshape(Xi, numel(idxPool), []);     % [cp x (Tw*N)]
    v  = var(double(Xi), 0, 2, 'omitnan');

    k  = max(1, round(frac * numel(idxPool)));
    [v_sorted, ord] = sort(v, 'ascend');
    sel = ord(1:k);

    out(b).car_channels = idxPool(sel);
    out(b).var_sel_thr  = v_sorted(min(k, numel(v_sorted)));
    out(b).n_sel        = numel(sel);
end

% ----- apply CAR per epoch with stim-channel exclusion --------------------
for b = 1:nBlk
    idxSet    = out(b).channels_set;
    idxCarBase= out(b).car_channels;
    if isempty(idxSet) || isempty(idxCarBase), continue; end
    for tr = 1:N
        if ~isempty(stim_excl{tr})
            idxCar = setdiff(idxCarBase, stim_excl{tr});
            if isempty(idxCar), idxCar = idxCarBase; end % fallback if all excluded
        else
            idxCar = idxCarBase;
        end
        mu = mean(X(idxCar,:,tr), 1, 'omitnan');   % [1 x T]
        X(idxSet,:,tr) = X(idxSet,:,tr) - mu;      % subtract broadcast
    end
end

% ----- write back & rank proxy -----
EEG.data = X;
EEG.ref  = 'aCAR';
r_after  = local_rank_proxy(X);
if r_after < r_before
    warning('aCAR: effective rank proxy decreased: %d -> %d', r_before, r_after);
end

% Print citation & docs pointer for aCAR in iEEGLAB
fprintf('\nPlease cite the following when performing aCAR with iEEGLAB:\n');
fprintf('Huang H., et al., (2024). CARLA: Adjusted common average referencing for cortico-cortical \n');
fprintf('evoked potential data. Journal of Neuroscience Methods, 407:110153. \n https://doi.org/10.1016/j.jneumeth.2024.110153\n');


% ======================= local helpers =======================
function r = local_rank_proxy(X3)
    X2 = reshape(double(X3), size(X3,1), []);  % C x (T*N)
    X2 = X2 - mean(X2,2,'omitnan');
    Cc = cov(X2');                             
    e  = eig(Cc);
    r  = sum(e > max(e)*1e-7);
end

function stim_excl = local_collect_stim_exclusions(EEG, chanLabels, Nexp)
    % Return cell{N} of channel indices to exclude per epoch
    stim_excl = repmat({[]}, 1, Nexp);
    map_one = @(tokens) local_tokens_to_idx(tokens, chanLabels);

    % Priority 1: EEG.ieeglab.opt.events (1 row per epoch)
    ev = [];
    if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'opt') && isfield(EEG.ieeglab.opt,'events')
        ev = EEG.ieeglab.opt.events;
        if istable(ev) && height(ev)~=Nexp
            ev = []; % size mismatch; ignore
        end
    end

    if ~isempty(ev)
        fields_try = {'electrodes_involved_onset','electrical_stimulation_site','stim_electrodes'};
        for i=1:Nexp
            tokens = [];
            for f = 1:numel(fields_try)
                if ismember(fields_try{f}, ev.Properties.VariableNames)
                    tokens = local_split_labels(ev{i,fields_try{f}}); 
                    if ~isempty(tokens), break; end
                end
            end
            if isempty(tokens)
                % Try stim_type + stim_num → like "LA" + "1" → "LA1"
                has_type = ismember('stim_type', ev.Properties.VariableNames);
                has_num  = ismember('stim_num',  ev.Properties.VariableNames);
                if has_type && has_num
                    t = ev.stim_type{i}; n = ev.stim_num(i);
                    try
                        tokens = cellstr(string(t) + string(n));
                    catch
                        tokens = {};
                    end
                end
            end
            stim_excl{i} = map_one(tokens);
        end
        return
    end

    % Priority 2: EEG.epoch.* (fields mirrored into epoch struct)
    if isfield(EEG,'epoch') && ~isempty(EEG.epoch) && numel(EEG.epoch)==Nexp
        fields_try = {'electrodes_involved_onset','electrical_stimulation_site','stim_electrodes'};
        for i=1:Nexp
            tokens = [];
            for f = 1:numel(fields_try)
                fn = fields_try{f};
                if isfield(EEG.epoch, fn)
                    tokens = local_split_labels(EEG.epoch(i).(fn));
                    if ~isempty(tokens), break; end
                end
            end
            % If nothing, try event fields joined into a string
            if isempty(tokens)
                if isfield(EEG.epoch,'stim_type') && isfield(EEG.epoch,'stim_num')
                    t = EEG.epoch(i).stim_type; n = EEG.epoch(i).stim_num;
                    try, tokens = cellstr(string(t) + string(n)); catch, end
                end
            end
            stim_excl{i} = map_one(tokens);
        end
    end
end

function tokens = local_split_labels(val)
    % Accept cell/char/string/numeric; split on common delimiters
    tokens = {};
    if isempty(val), return; end
    if iscell(val), val = string(val); end
    if iscategorical(val), val = string(val); end
    if isnumeric(val)
        tokens = cellstr(string(val));
        return
    end
    if isstring(val)
        s = strjoin(val, ' ');
        parts = regexp(char(s), '[,;+\-\/\|\s]+', 'split');
        tokens = parts(~cellfun(@isempty,parts));
    elseif ischar(val)
        parts = regexp(val, '[,;+\-\/\|\s]+', 'split');
        tokens = parts(~cellfun(@isempty,parts));
    end
end

function idx = local_tokens_to_idx(tokens, chanLabels)
    idx = [];
    if isempty(tokens), return; end
    tok = string(tokens(:));
    [tf,loc] = ismember(upper(tok), upper(chanLabels)); % exact matches
    idx = loc(tf);
    if isempty(idx)
        tok2 = regexprep(tok, '[^\w]', '');             % strip non-word chars
        [tf2,loc2] = ismember(upper(tok2), upper(chanLabels));
        idx = loc2(tf2);
    end
    idx = unique(idx(idx>0));
end
end
