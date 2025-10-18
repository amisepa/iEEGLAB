function EEG = get_elec_coor(EEG, elecs)
%GET_ELEC_COOR  Copy coordinates from an electrodes TSV/table into EEG.chanlocs.
%
% Matching strategy per channel:
%   1) Exact (case-insensitive) match between EEG.chanlocs(i).labels and elecs{:,ilab}
%   2) Fallback: extract integer from both labels (e.g., '1' <-> 'iEEG1')
%      and match by that number (first unique match wins)
%
% Notes:
% - The first column of the 'elecs' table is treated as the electrode label column.
% - X/Y/Z are assigned only if present and finite. NaNs are skipped.
% - If a fallback numeric match is used, the matched TSV name is stored in
%   EEG.chanlocs(i).matched_elec_label for traceability.

assert(istable(elecs), 'elecs must be a table');

% Column indices (case-insensitive)
vars = elecs.Properties.VariableNames;
ilab = 1;  % first column holds the channel/electrode labels
ix   = find(strcmpi(vars,'x'), 1);
iy   = find(strcmpi(vars,'y'), 1);
iz   = find(strcmpi(vars,'z'), 1);
idn  = find(strcmpi(vars,'Destrieux_label'), 1);        % number/index (optional)
idl  = find(strcmpi(vars,'Destrieux_label_text'), 1);   % text label (optional)

% Normalize electrode names column to string for consistent handling
elec_names = string(elecs{:, ilab});

% Precompute lowercase versions for exact/case-insensitive matching
elec_names_lc = lower(strtrim(elec_names));

% Precompute numeric ids extracted from electode names (e.g., "iEEG114" -> 114)
elec_nums = extract_first_integer(elec_names);

missing = {};
for iChan = 1:EEG.nbchan
    lab  = string(EEG.chanlocs(iChan).labels);
    labL = lower(strtrim(lab));

    r = [];  % row index into elecs

    % 1) Try exact/case-insensitive text match
    hit = find(elec_names_lc == labL, 1, 'first');
    if ~isempty(hit)
        r = hit;
    else
        % 2) Fallback: numeric match if both sides expose an integer
        ch_num = extract_first_integer(lab);
        if ~isnan(ch_num)
            % Find first occurrence of the same number among elec labels
            hit_num = find(elec_nums == ch_num, 1, 'first');
            if ~isempty(hit_num)
                r = hit_num;
            end
        end
    end

    if ~isempty(r)
        % Assign coordinates if available and finite
        if ~isempty(ix), EEG.chanlocs(iChan).X = safe_num(elecs{r, ix}); end
        if ~isempty(iy), EEG.chanlocs(iChan).Y = safe_num(elecs{r, iy}); end
        if ~isempty(iz), EEG.chanlocs(iChan).Z = safe_num(elecs{r, iz}); end

        % Optional extras if present
        if ~isempty(idn) && ~ismissing(elecs{r, idn})
            EEG.chanlocs(iChan).destrieux_num = elecs{r, idn};
        end
        if ~isempty(idl) && ~ismissing(elecs{r, idl})
            val = elecs{r, idl};
            EEG.chanlocs(iChan).destrieux_label = char(string(val));
        end

        % Record which TSV label we matched to (useful when numeric fallback used)
        EEG.chanlocs(iChan).matched_elec_label = char(elec_names(r));
    else
        missing{end+1} = char(lab); %#ok<AGROW>
    end
end

if ~isempty(missing)
    warning('No electrode coordinates found for %d channel(s): %s', ...
        numel(missing), strjoin(missing, ', '));
end
end

% ---------- helpers ----------
function v = safe_num(x)
% Convert to double if finite; otherwise leave empty (avoid writing NaN).
    if isnumeric(x) && isscalar(x) && isfinite(x)
        v = double(x);
    else
        v = [];
    end
end

function n = extract_first_integer(s)
% Return the first integer found in a string/char/cellstr; NaN if none.
    if iscell(s), s = string(s); end
    if ~isstring(s), s = string(s); end
    n = NaN(size(s));
    for k = 1:numel(s)
        tk = char(s(k));
        % Grab the first run of digits
        m = regexp(tk, '\d+', 'match', 'once');
        if ~isempty(m)
            n(k) = str2double(m);
        end
    end
end
