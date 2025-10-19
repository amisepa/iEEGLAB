function EEG = get_elec_coor(EEG, elecs)
%GET_ELEC_COOR  Copy X/Y/Z from table into EEG.chanlocs with robust matching.
%
% Attempts:
%   1) Direct case-insensitive exact label match.
%   2) If 0% matched: normalize BOTH EEG and TSV labels by removing "EEG"/"eeg"
%      and tidying separators (without mutating visible labels), then retry.
%   3) If still 0% matched: overwrite EEG.chanlocs(i).labels from TSV order, then retry.
%
% Also:
%   - Converts X/Y/Z from strings/cells to doubles safely.
%   - Prints brief sanity previews (label changes, coordinates).

assert(istable(elecs), 'elecs must be a table');

% ---------- Column indices (case-insensitive) ----------
ilab = 1;
vars = elecs.Properties.VariableNames;
ix = find(strcmpi(vars,'x'), 1);
iy = find(strcmpi(vars,'y'), 1);
iz = find(strcmpi(vars,'z'), 1);
assert(~isempty(ix) && ~isempty(iy) && ~isempty(iz), ...
    'The electrodes table must have X, Y, Z columns (case-insensitive).');

% TSV label lookup
tsv_labels    = string(elecs{:, ilab});
tsv_labels_lc = lower(strtrim(tsv_labels));

% ---------- Attempt 1: direct label match ----------
[EEG, matchedMask] = copy_xyz_if_matched(EEG, elecs, tsv_labels, tsv_labels_lc, ix, iy, iz);
report_match('[get_elec_coor] Attempt 1 (as-is)', EEG, matchedMask);

% ---------- Attempt 2: normalized comparison (remove "EEG"/tidy) w/o relabel ----------
if sum(matchedMask) == 0
    eeg_labels        = string({EEG.chanlocs.labels}.');
    eeg_labels_norm   = normalize_label(eeg_labels);        % remove "EEG", tidy
    tsv_labels_norm   = normalize_label(tsv_labels);        % same transform on TSV
    tsv_labels_norm_l = lower(strtrim(tsv_labels_norm));

    [EEG, matchedMask] = copy_xyz_if_matched_norm( ...
        EEG, elecs, tsv_labels, tsv_labels_norm_l, eeg_labels_norm, ix, iy, iz);

    report_match('[get_elec_coor] Attempt 2 (normalized "no-EEG")', EEG, matchedMask);
end

% ---------- Attempt 3: overwrite labels from TSV order and retry ----------
if sum(matchedMask) == 0
    oldLabels = string({EEG.chanlocs.labels}.');
    m = min(EEG.nbchan, height(elecs));
    if m < EEG.nbchan
        warning('TSV has fewer labels (%d) than EEG channels (%d). Only first %d channels will be relabeled.', ...
            m, EEG.nbchan, m);
    end

    % Overwrite labels in order
    for i = 1:m
        EEG.chanlocs(i).labels = char(tsv_labels(i));
    end

    % Sanity print: first 20 (or all if <20)
    k = min(20, m);
    fprintf('[get_elec_coor] Attempt 3 relabel: first %d channel(s) from TSV order (old -> new):\n', k);
    for i = 1:k
        fprintf('  [%3d] "%s" -> "%s"\n', i, oldLabels(i), string(EEG.chanlocs(i).labels));
    end
    if m > k
        fprintf('  ... and %d more relabeled.\n', m - k);
    end

    % Re-attempt coordinate copy after relabel
    [EEG, matchedMask] = copy_xyz_if_matched(EEG, elecs, tsv_labels, tsv_labels_lc, ix, iy, iz);
    report_match('[get_elec_coor] After Attempt 3 (TSV relabel)', EEG, matchedMask);
end

% ---------- Coordinate sanity preview ----------
preview_n = min(5, EEG.nbchan);
fprintf('[get_elec_coor] First %d coordinate rows after assignment:\n', preview_n);
for ii = 1:preview_n
    c = EEG.chanlocs(ii);
    fprintf('  [%3d] %-12s  X=%s  Y=%s  Z=%s\n', ii, c.labels, ...
        num2str_safe(c, 'X'), num2str_safe(c, 'Y'), num2str_safe(c, 'Z'));
end

end

% ---- helpers ----
function [EEG, matchedMask] = copy_xyz_if_matched(EEG, elecs, tsv_labels, tsv_labels_lc, ix, iy, iz)
matchedMask = false(EEG.nbchan,1);
for i = 1:EEG.nbchan
    key = lower(strtrim(string(EEG.chanlocs(i).labels)));
    hit = find(tsv_labels_lc == key, 1, 'first');
    if ~isempty(hit)
        EEG.chanlocs(i).matched_elec_label = char(tsv_labels(hit));
        xi = toNumericScalar(elecs{hit, ix});
        yi = toNumericScalar(elecs{hit, iy});
        zi = toNumericScalar(elecs{hit, iz});
        if ~isnan(xi), EEG.chanlocs(i).X = xi; end
        if ~isnan(yi), EEG.chanlocs(i).Y = yi; end
        if ~isnan(zi), EEG.chanlocs(i).Z = zi; end
        matchedMask(i) = true;
    end
end
end

function [EEG, matchedMask] = copy_xyz_if_matched_norm(EEG, elecs, tsv_labels_orig, tsv_labels_norm_l, eeg_labels_norm, ix, iy, iz)
% Compare normalized EEG labels to normalized TSV labels, but assign XYZ
% from the original TSV row (tsv_labels_orig). Does NOT change visible labels.
matchedMask = false(EEG.nbchan,1);
for i = 1:EEG.nbchan
    key_norm = lower(strtrim(string(eeg_labels_norm(i))));
    hit = find(tsv_labels_norm_l == key_norm, 1, 'first');
    if ~isempty(hit)
        EEG.chanlocs(i).matched_elec_label = char(tsv_labels_orig(hit));
        xi = toNumericScalar(elecs{hit, ix});
        yi = toNumericScalar(elecs{hit, iy});
        zi = toNumericScalar(elecs{hit, iz});
        if ~isnan(xi), EEG.chanlocs(i).X = xi; end
        if ~isnan(yi), EEG.chanlocs(i).Y = yi; end
        if ~isnan(zi), EEG.chanlocs(i).Z = zi; end
        matchedMask(i) = true;
    end
end
end

function lab = normalize_label(lab)
% Remove any "EEG"/"eeg" substring and tidy separators; keeps original case of others.
    lab = regexprep(string(lab), '(?i)eeg', '');   % strip "EEG" case-insensitive
    lab = strip_separators(lab);                   % trim ends, collapse runs
end

function out = strip_separators(lbl)
% Trim leading/trailing separators/spaces and collapse repeats to a single underscore.
    out = regexprep(string(lbl), '(^[\s_\-]+|[\s_\-]+$)', '');
    out = regexprep(out, '[\s_\-]{2,}', '_');
    out = strtrim(out);
end

function report_match(tag, EEG, matchedMask)
nTotal     = EEG.nbchan;
nMatched   = sum(matchedMask);
nUnmatched = nTotal - nMatched;
pctMatched = 100 * nMatched / max(nTotal,1);
fprintf('%s: Matched %d / %d (%.1f%%). Unmatched: %d (%.1f%%).\n', ...
    tag, nMatched, nTotal, pctMatched, nUnmatched, 100 - pctMatched);
end

function val = toNumericScalar(x)
% Returns a double scalar or NaN. Handles cell, string, char, missing.
    if iscell(x), x = x{1}; end
    if isstring(x), x = char(x); end
    if ischar(x)
        x = strrep(x, ',', '.');    % allow comma decimals
        x = str2double(strtrim(x));
    end
    if isnumeric(x) && isscalar(x)
        val = double(x);
    else
        val = NaN;
    end
end

function s = num2str_safe(c, fld)
    if isfield(c, fld) && ~isempty(c.(fld))
        s = num2str(c.(fld));
    else
        s = '[]';
    end
end
