function EEG = get_elec_coor(EEG, elecs)

assert(istable(elecs), 'elecs must be a table');

% Column indices (case-insensitive); OK if some are missing
vars = elecs.Properties.VariableNames;
ilab = 1;  % first column holds the channel/electrode labels
ix   = find(strcmpi(vars,'x'), 1);
iy   = find(strcmpi(vars,'y'), 1);
iz   = find(strcmpi(vars,'z'), 1);
idn  = find(strcmpi(vars,'Destrieux_label'), 1);        % number/index
idl  = find(strcmpi(vars,'Destrieux_label_text'), 1);   % text label

missing = {};
for iChan = 1:EEG.nbchan
    lab = EEG.chanlocs(iChan).labels;
    match = strcmpi(elecs{:, ilab}, lab);

    if any(match)
        r = find(match, 1, 'first');

        if ~isempty(ix) && ~ismissing(elecs{r, ix})
            EEG.chanlocs(iChan).X = double(elecs{r, ix});
        end
        if ~isempty(iy) && ~ismissing(elecs{r, iy})
            EEG.chanlocs(iChan).Y = double(elecs{r, iy});
        end
        if ~isempty(iz) && ~ismissing(elecs{r, iz})
            EEG.chanlocs(iChan).Z = double(elecs{r, iz});
        end
        if ~isempty(idn) && ~ismissing(elecs{r, idn})
            EEG.chanlocs(iChan).destrieux_num = elecs{r, idn};
        end
        if ~isempty(idl) && ~ismissing(elecs{r, idl})
            val = elecs{r, idl};
            EEG.chanlocs(iChan).destrieux_label = char(string(val));
        end
    else
        missing{end+1} = char(lab); %#ok<AGROW>
    end
end

if ~isempty(missing)
    warning('No electrode coordinates found for %d channel(s): %s', ...
        numel(missing), strjoin(missing, ', '));
end

