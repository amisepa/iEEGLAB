function [opt, wasCanceled] = ieeglab_gui_load2(opt)
% iEEGLAB GUI (channels, minimal event-based loading)
%
% Required:
%   opt.elec_labels  - cellstr of channel labels
%
% Optional (preselects only):
%   opt.chan_idx     - previously selected channel indices (preselect)
%   opt.events       - table; if present enables event-based loading UI
%   opt.event_field  - (optional) preselect events column name
%   opt.event_values - (optional) preselect values (cell of char) for that column
%
% Outputs set by GUI:
%   opt.chan_idx, opt.chan_list
%   opt.load_by_events       - logical
%   opt.event_field          - char (when load_by_events==true, else '')
%   opt.event_values         - {'All'} or list of values (cell of char) (when load_by_events)
%   opt.epoch_window_sec     - [pre post] seconds (default [-1 1.5]) (when load_by_events)
%   wasCanceled  logical, true if user canceled

wasCanceled = false;      % default

assert(isfield(opt,'elec_labels') && ~isempty(opt.elec_labels), ...
    'ieeglab_gui_load2: opt.elec_labels must be a non-empty cell array.');
labels = opt.elec_labels(:);
nch    = numel(labels);

% ---------- Channels preselect ----------
if ~isfield(opt,'chan_idx') || isempty(opt.chan_idx)
    preCh = 1;
else
    sel = opt.chan_idx(:)'; sel = sel(sel>=1 & sel<=nch);
    preCh = iff(isempty(sel) || numel(sel)==nch, 1, sel+1);
end
labels_disp = [{'All channels'}; labels];

% ---------- Events presence & lists ----------
haveEvents = isfield(opt,'events') && istable(opt.events) && ~isempty(opt.events);

fields_disp = {};
all_value_lists = {};
preFieldIdx = 1;

if haveEvents
    vars = opt.events.Properties.VariableNames;
    cand = setdiff(vars, {'onset'}, 'stable');  % all event fields except onset
    for k = 1:numel(cand)
        col  = cand{k};
        vals = getUniques(opt.events, col);                   % string vector
        fields_disp{end+1}     = col;                         %#ok<AGROW>
        all_value_lists{end+1} = [{'All'}; cellstr(string(vals(:)))]; %#ok<AGROW>
    end
    if isempty(fields_disp)
        fields_disp = {'<no event fields>'};
        all_value_lists = {{'All'}};
    end
    if isfield(opt,'event_field') && ~isempty(opt.event_field)
        fidx = find(strcmp(fields_disp, opt.event_field), 1);
        if ~isempty(fidx), preFieldIdx = fidx; end
    else
        def = find(strcmp(fields_disp,'electrical_stimulation_site'),1);
        if ~isempty(def), preFieldIdx = def; end
    end
else
    fields_disp     = {'<no events table>'};
    all_value_lists = {{'All'}};
    preFieldIdx     = 1;
end

% Values list for current field (used when load_by_events is on)
list_values = all_value_lists{preFieldIdx};

% Preselect values (from opt.event_values) for the listbox (if later enabled)
preValIdx = 1; % "All"
if isfield(opt,'event_values') && ~isempty(opt.event_values)
    want = cellstr(opt.event_values(:));
    idxs = [];
    for w = 1:numel(want)
        i = find(strcmp(list_values, want{w}), 1);
        if ~isempty(i), idxs(end+1) = i; end %#ok<AGROW>
    end
    if ~isempty(idxs), preValIdx = idxs; end
end

% ---------- Callbacks ----------
% Toggle event block on/off (by tags)
cb_load_by_ev = [
 'val=get(gcbo,''value'');' ...
 'if val, onoff=''on''; else, onoff=''off''; end;' ...
 'set(findobj(gcbf,''tag'',''ev_field_label''),''enable'',onoff);' ...
 'set(findobj(gcbf,''tag'',''event_field''),''enable'',onoff);' ...
 'set(findobj(gcbf,''tag'',''ev_values_label''),''enable'',onoff);' ...
 'set(findobj(gcbf,''tag'',''event_val_list''),''enable'',onoff);' ...
 'set(findobj(gcbf,''tag'',''ev_epoch_label''),''enable'',onoff);' ...
 'set(findobj(gcbf,''tag'',''epoch_win''),''enable'',onoff);' ...
];

% Update values list when event-field popup changes (vallists stored in popup UserData)
cb_field = [
 'idx = get(gcbo,''value'');' ...
 'vallists = get(gcbo,''userdata'');' ...
 'lst = vallists{idx};' ...
 'hlist = findobj(gcbf,''tag'',''event_val_list'');' ...
 'set(hlist,''string'',lst,''value'',1);' ...
];

% ---------- Geometry (7 rows) ----------
% 1) text, 2) listbox, 3) "Loading options", 4) checkbox,
% 5) field label+popup, 6) values label+listbox, 7) epoch label+edit
uigeom = {
    [1]      % row 1
    [1]      % row 2
    [1]      % row 3
    [1]      % row 4
    [1 1]    % row 5
    [1 1]    % row 6
    [1 1]    % row 7
};
geomvert = [1 6 1 1 1 2 1];  % heights per row (must match rows in uigeom)

% ---------- UI list (10 components total, must match uigeom columns) ----------
uilist = {
    % Row 1: Channels title
    {'style' 'text'    'string' 'Select channel(s):' 'fontweight' 'bold'}
    % Row 2: Channels listbox
    {'style' 'listbox' 'tag' 'sel_list' 'string' labels_disp ...
     'value' preCh 'min' 0 'max' 2 'ListboxTop' 1}

    % Row 3: Loading options header
    {'style' 'text' 'string' 'Loading options' 'fontweight' 'bold'}

    % Row 4: Load by events (renamed)
    {'style' 'checkbox' 'tag' 'load_by_events' ...
     'string' 'Load data by events' ...
     'value' 0 'callback' cb_load_by_ev 'enable' iff(haveEvents,'on','off')}

    % Row 5: Event field picker (label + popup)
    {'style' 'text' 'tag' 'ev_field_label' 'string' 'Event field:' 'enable' 'off'}
    {'style' 'popupmenu' 'tag' 'event_field' 'string' fields_disp ...
     'value' preFieldIdx 'callback' cb_field ...
     'userdata' all_value_lists 'enable' 'off'}

    % Row 6: Values for that field (label + listbox)
    {'style' 'text' 'tag' 'ev_values_label' 'string' 'Values to load:' 'enable' 'off'}
    {'style' 'listbox' 'tag' 'event_val_list' 'string' list_values ...
     'value' preValIdx 'min' 0 'max' 2 'ListboxTop' 1 'enable' 'off'}

    % Row 7: Epoch window (label + edit)
    {'style' 'text' 'tag' 'ev_epoch_label' 'string' 'Epoch window (s) [pre post]:' 'enable' 'off'}
    {'style' 'edit' 'tag' 'epoch_win' 'string' '-0.5 1' ...
     'horizontalalignment' 'left' 'enable' 'off'}
};

[res,~,~,outstruct] = inputgui('geometry', uigeom, ...
                             'geomvert', geomvert, ...
                             'uilist', uilist, ...
                             'title','iEEGLAB: load options', ...
                             'minwidth', 600);
% Clean abort on Cancel or window close
if isempty(res) || isempty(outstruct)
    fprintf('iEEGLAB load dialog canceled\n');
    wasCanceled = true;   % mark canceled
    opt = [];
    return;
end

% ---------- Channels mapping ----------
idxCh_disp = outstruct.sel_list(:)';
if any(idxCh_disp==1) || isempty(idxCh_disp)
    opt.chan_idx  = 1:nch;
    opt.chan_list = labels;
else
    idx           = idxCh_disp - 1;
    opt.chan_idx  = idx;
    opt.chan_list = labels(idx);
end

% ---------- Outputs for event-based loading ----------
opt.load_by_events   = isfield(outstruct,'load_by_events') && outstruct.load_by_events==1;
if opt.load_by_events && haveEvents
    opt.event_field      = '';
    opt.event_values     = {'All'};
    opt.epoch_window_sec = [-1 1.5];

    % event field index (with fallback)
    if isfield(outstruct,'event_field') && ~isempty(outstruct.event_field)
        ef_idx = outstruct.event_field;
    else
        ef_idx = preFieldIdx;
    end
    ef_idx = max(1, min(ef_idx, numel(fields_disp)));
    opt.event_field = fields_disp{ef_idx};

    % values list (with fallback)
    vals_for_field = all_value_lists{ef_idx};
    if isfield(outstruct,'event_val_list') && ~isempty(outstruct.event_val_list)
        opt.event_values = mapOut(outstruct.event_val_list, vals_for_field);
    else
        opt.event_values = {'All'};
    end

    % epoch window
    if isfield(outstruct,'epoch_win') && ~isempty(outstruct.epoch_win)
        ew_raw  = strtrim(outstruct.epoch_win);
        ew_nums = regexp(ew_raw, '([-+]?\d*\.?\d+)', 'match');
        if numel(ew_nums) >= 2
            pre  = str2double(ew_nums{1});
            post = str2double(ew_nums{2});
            if ~isnan(pre) && ~isnan(post), opt.epoch_window_sec = [pre post]; end
        end
    end
else
    % % Not loading by events
    % opt.event_field  = '';
    % opt.event_values = {'All'};
end

end % === main ===

% -------------------- helpers --------------------
function out = mapOut(idx_list, list_disp)
    % Return cell array of char (row cell). 'All' => {'All'}
    if isempty(idx_list) || any(idx_list==1)
        out = {'All'};
    else
        out = list_disp(idx_list);   % list_disp is cellstr
        out = out(:)';               % row cell
    end
end

function vals = getUniques(ev, colname)
    if ismember(colname, ev.Properties.VariableNames)
        x = ev.(colname);
        if iscell(x), x = string(x); end
        vals = unique(x(~ismissing(x)));
    else
        vals = string.empty(0,1);
    end
end

function y = iff(cond, a, b)
    if cond, y = a; else, y = b; end
end
