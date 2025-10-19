function EEG = get_elec_coor(EEG, elecs)
%GET_ELEC_COOR  Exact label match only: copy X/Y/Z from table into EEG.chanlocs.
%
% Behavior
%   - For each EEG channel i:
%       * Find case-insensitive exact match to the first column of 'elecs'.
%       * If found, copy X/Y/Z (if finite) into EEG.chanlocs(i).
%       * If not found, leave EEG.chanlocs(i) unchanged and report later.
%
% Assumptions
%   - 'elecs' is a table; first column is the label column.
%   - Coordinate columns are named 'X','Y','Z' (any case).

assert(istable(elecs), 'elecs must be a table');

% Column indices
ilab = 1;
vars = elecs.Properties.VariableNames;
ix = find(strcmpi(vars,'x'), 1);
iy = find(strcmpi(vars,'y'), 1);
iz = find(strcmpi(vars,'z'), 1);
assert(~isempty(ix) && ~isempty(iy) && ~isempty(iz), ...
    'The electrodes table must have X, Y, Z columns (case-insensitive).');

% Prepare label lookup (lowercased / trimmed)
tsv_labels = string(elecs{:, ilab});
tsv_labels_lc = lower(strtrim(tsv_labels));

missing = strings(0,1);
matchedMask = false(EEG.nbchan,1);

for i = 1:EEG.nbchan
    lab = string(EEG.chanlocs(i).labels);
    key = lower(strtrim(lab));

    % exact, case-insensitive match (take the first if duplicated)
    hit = find(tsv_labels_lc == key, 1, 'first');

    if ~isempty(hit)
        EEG.chanlocs(i).matched_elec_label = char(tsv_labels(hit));
        xi = elecs{hit, ix}; yi = elecs{hit, iy}; zi = elecs{hit, iz};
        if isfinite_scalar(xi), EEG.chanlocs(i).X = double(xi); end
        if isfinite_scalar(yi), EEG.chanlocs(i).Y = double(yi); end
        if isfinite_scalar(zi), EEG.chanlocs(i).Z = double(zi); end
        matchedMask(i) = true;
    else
        missing(end+1,1) = lab; %#ok<AGROW>
    end
end

% Summary reporting
nTotal    = EEG.nbchan;
nMatched  = sum(matchedMask);
nUnmatched = nTotal - nMatched;
pctMatched = 100 * nMatched / max(nTotal,1);

fprintf('[get_elec_coor] Matched: %d / %d (%.1f%%). Unmatched: %d (%.1f%%).\n', ...
    nMatched, nTotal, pctMatched, nUnmatched, 100 - pctMatched);

if nUnmatched > 0
    warning('No XYZ match for %d channel(s): %s', ...
        nUnmatched, strjoin(cellstr(missing.'), ', '));
end
end

% ---- helpers ----
function tf = isfinite_scalar(x)
    tf = isnumeric(x) && isscalar(x) && isfinite(x);
end
