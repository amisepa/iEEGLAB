function [EEG, wasCanceled] = ieeglab_gui_preprocess(EEG)

wasCanceled = false;

% dynamic defaults
ds_should_enable = EEG.srate > 512;
ds_default_rate  = 512;

% defaults
choices = struct();
choices.event_filters     = struct();
choices.remove_rare_cond  = true;
choices.min_trials        = 5;
choices.remove_no_coords  = true;

choices.apply_ds          = ds_should_enable;
choices.apply_highpass    = true;
choices.apply_notch       = true;
choices.apply_lowpass     = false;
choices.apply_epoch       = true;

if ds_should_enable
    choices.downsample    = ds_default_rate;
else
    choices.downsample    = max(1, EEG.srate);
end
choices.highpass          = 0.1;
choices.highpass_type     = 1;      % 1=noncausal zero-phase, 2=causal minimum-phase
choices.notch             = 60;
choices.lowpass           = [];
choices.epoch_window      = [-500 900]; % ms

choices.apply_acar        = true;
choices.acar_fraction     = 0.20;
choices.acar_timewin      = [15 500];   % ms

choices.apply_baseline    = true;
choices.baseline_method   = 1;      % 1=median, 2=mean, 3=trimmed mean, 4=1/f
choices.baseline_period   = [-500 -50]; % ms
choices.baseline_mode     = 1;      % 1=subtract, 2=divide

% events availability
haveEv = isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'opt') && ...
         isfield(EEG.ieeglab.opt,'events') && istable(EEG.ieeglab.opt.events) && ...
         ~isempty(EEG.ieeglab.opt.events);

if haveEv
    evT = EEG.ieeglab.opt.events;
    ignoreCols = {'onset','duration','electrodes_involved_onset','electrodes_involved_offset','sample_start'};
    allCols    = evT.Properties.VariableNames;
    showCols   = setdiff(allCols, ignoreCols, 'stable');
    mustInclude = intersect({'stim_type','stim_num'}, allCols, 'stable');
    showCols   = unique([showCols(:); mustInclude(:)], 'stable');
else
    evT = table(); showCols = {};
end

% callbacks
if haveEv && ~isempty(showCols)
    cb_ev = @(h,~) local_open_event_selector(h, evT, showCols);
    btn_en = 'on';
else
    cb_ev = @(~,~) warndlg('No events table available.','iEEGLAB');
    btn_en = 'off';
end

% cb_ds  = ['val=get(gcbo,''value''); onoff=tern(val);' ...
%           'set(findobj(gcbf,''tag'',''lbl_ds_rate''),''enable'',onoff);' ...
%           'set(findobj(gcbf,''tag'',''ds_rate''),''enable'',onoff);'];
% cb_hp  = ['val=get(gcbo,''value''); onoff=tern(val);' ...
%           'set(findobj(gcbf,''tag'',''lbl_highpass''),''enable'',onoff);' ...
%           'set(findobj(gcbf,''tag'',''highpass''),''enable'',onoff);' ...
%           'set(findobj(gcbf,''tag'',''lbl_highpass_type''),''enable'',onoff);' ...
%           'set(findobj(gcbf,''tag'',''highpass_type''),''enable'',onoff);'];
% cb_no  = ['val=get(gcbo,''value''); onoff=tern(val);' ...
%           'set(findobj(gcbf,''tag'',''lbl_notch''),''enable'',onoff);' ...
%           'set(findobj(gcbf,''tag'',''notch''),''enable'',onoff);'];
% cb_lp  = ['val=get(gcbo,''value''); onoff=tern(val);' ...
%           'set(findobj(gcbf,''tag'',''lbl_lowpass''),''enable'',onoff);' ...
%           'set(findobj(gcbf,''tag'',''lowpass''),''enable'',onoff);'];
% cb_seg = ['val=get(gcbo,''value''); onoff=tern(val);' ...
%           'set(findobj(gcbf,''tag'',''lbl_epoch''),''enable'',onoff);' ...
%           'set(findobj(gcbf,''tag'',''epoch_window''),''enable'',onoff);'];
% cb_acar= ['val=get(gcbo,''value''); onoff=tern(val);' ...
%           'set(findobj(gcbf,''tag'',''acar_fraction''),''enable'',onoff);' ...
%           'set(findobj(gcbf,''tag'',''acar_timewin''),''enable'',onoff);' ...
%           'set(findobj(gcbf,''tag'',''lbl_acar_fraction''),''enable'',onoff);' ...
%           'set(findobj(gcbf,''tag'',''lbl_acar_timewin''),''enable'',onoff);'];
% cb_bl  = ['val=get(gcbo,''value''); onoff=tern(val);' ...
%           'set(findobj(gcbf,''tag'',''lbl_bl_method''),''enable'',onoff);' ...
%           'set(findobj(gcbf,''tag'',''bl_method''),''enable'',onoff);' ...
%           'set(findobj(gcbf,''tag'',''lbl_bl_period''),''enable'',onoff);' ...
%           'set(findobj(gcbf,''tag'',''bl_period''),''enable'',onoff);' ...
%           'set(findobj(gcbf,''tag'',''lbl_bl_mode''),''enable'',onoff);' ...
%           'set(findobj(gcbf,''tag'',''bl_mode''),''enable'',onoff);'];
cb_ds  = local_cb_toggle({'lbl_ds_rate','ds_rate'});
cb_hp  = local_cb_toggle({'lbl_highpass','highpass','lbl_highpass_type','highpass_type'});
cb_no  = local_cb_toggle({'lbl_notch','notch'});
cb_lp  = local_cb_toggle({'lbl_lowpass','lowpass'});
cb_seg = local_cb_toggle({'lbl_epoch','epoch_window'});
cb_acar= local_cb_toggle({'acar_fraction','acar_timewin','lbl_acar_fraction','lbl_acar_timewin'});
cb_bl  = local_cb_toggle({'lbl_bl_method','bl_method','lbl_bl_period','bl_period','lbl_bl_mode','bl_mode'});

onoff_ds  = iff(choices.apply_ds,'on','off');
onoff_hp  = iff(choices.apply_highpass,'on','off');
onoff_no  = iff(choices.apply_notch,'on','off');
onoff_lp  = iff(choices.apply_lowpass,'on','off');
onoff_seg = iff(choices.apply_epoch,'on','off');
onoff_car = iff(choices.apply_acar,'on','off');
onoff_bl  = iff(choices.apply_baseline,'on','off');

%% geometry (narrower; "remove conditions" enable + indented value)
uigeom = {
    1
    [0.70 0.30]
    [0.70 0.30]
    [0.06 0.64 0.30]
    0.01
    1
    [0.70 0.30]
    1

    [0.70 0.30]
    [0.06 0.64 0.30]

    [0.70 0.30]
    [0.06 0.64 0.30]
    [0.06 0.64 0.30]

    [0.70 0.30]
    [0.06 0.64 0.30]

    [0.70 0.30]
    [0.06 0.64 0.30]

    [0.70 0.30]
    [0.06 0.64 0.30]

    [0.70 0.30]
    [0.06 0.64 0.30]
    [0.06 0.64 0.30]

    [0.70 0.30]
    [0.06 0.64 0.30]
    [0.06 0.64 0.30]
    [0.06 0.64 0.30]
};

%% build uilist with safe appends (no vertcat)
uilist = {};

% helper to append one control
    function add(ctrl)
        uilist{end+1} = ctrl;
    end

% Events
add({'style' 'text' 'string' 'Events' 'fontweight' 'bold' 'horizontalalignment' 'left'});
add({'style' 'text' 'string' 'Select events of interest:' 'horizontalalignment' 'left'});
add({'style' 'pushbutton' 'string' 'Open selector…' 'callback' cb_ev 'enable' btn_en});

add({'style' 'text' 'string' 'Remove rare conditions:' 'horizontalalignment' 'left'});
add({'style' 'checkbox' 'tag' 'remove_rare_cond' 'value' choices.remove_rare_cond 'string' 'Enable'});

add({'style' 'text' 'string' ''});
add({'style' 'text' 'string' 'Min trials:' 'horizontalalignment' 'left'});
add({'style' 'edit' 'tag' 'min_trials' 'string' num2str(choices.min_trials)});

% hidden store (tiny height row above Electrodes)
add({'style' 'edit' 'tag' 'evsel_json' 'string' '' 'visible' 'off'});

% Electrodes
add({'style' 'text' 'string' 'Electrodes' 'fontweight' 'bold' 'horizontalalignment' 'left'});
add({'style' 'text' 'string' 'Remove electrodes with no XYZ coordinates?' 'horizontalalignment' 'left'});
add({'style' 'checkbox' 'tag' 'rm_no_coords' 'value' choices.remove_no_coords 'string' ''});

% Signal Processing
add({'style' 'text' 'string' 'Signal Processing' 'fontweight' 'bold' 'horizontalalignment' 'left'});

% Downsample
add({'style' 'text' 'string' 'Downsample:' 'horizontalalignment' 'left'});
add({'style' 'checkbox' 'tag' 'apply_ds' 'value' choices.apply_ds 'string' 'Enable' 'callback' cb_ds});
add({'style' 'text' 'string' ''});
add({'style' 'text' 'tag' 'lbl_ds_rate' 'string' 'Target rate (Hz):' 'horizontalalignment' 'left' 'enable' onoff_ds});
add({'style' 'edit' 'tag' 'ds_rate' 'string' num2str(choices.downsample) 'enable' onoff_ds});

% High-pass
add({'style' 'text' 'string' 'High-pass filter:' 'horizontalalignment' 'left'});
add({'style' 'checkbox' 'tag' 'apply_highpass' 'value' choices.apply_highpass 'string' 'Enable' 'callback' cb_hp});
add({'style' 'text' 'string' ''});
add({'style' 'text' 'tag' 'lbl_highpass' 'string' 'Cutoff (Hz):' 'horizontalalignment' 'left' 'enable' onoff_hp});
add({'style' 'edit' 'tag' 'highpass' 'string' num2str(choices.highpass) 'enable' onoff_hp});
add({'style' 'text' 'string' ''});
add({'style' 'text' 'tag' 'lbl_highpass_type' 'string' 'Type:' 'horizontalalignment' 'left' 'enable' onoff_hp});
add({'style' 'popupmenu' 'tag' 'highpass_type' 'string' {'Noncausal zero-phase (default)' 'Causal minimum-phase'} ...
     'value' choices.highpass_type 'enable' onoff_hp});

% Notch
add({'style' 'text' 'string' 'Notch filter:' 'horizontalalignment' 'left'});
add({'style' 'checkbox' 'tag' 'apply_notch' 'value' choices.apply_notch 'string' 'Enable' 'callback' cb_no});
add({'style' 'text' 'string' ''});
add({'style' 'text' 'tag' 'lbl_notch' 'string' 'Center (Hz):' 'horizontalalignment' 'left' 'enable' onoff_no});
add({'style' 'edit' 'tag' 'notch' 'string' num2str(choices.notch) 'enable' onoff_no});

% Low-pass
add({'style' 'text' 'string' 'Low-pass filter:' 'horizontalalignment' 'left'});
add({'style' 'checkbox' 'tag' 'apply_lowpass' 'value' choices.apply_lowpass 'string' 'Enable' 'callback' cb_lp});
add({'style' 'text' 'string' ''});
add({'style' 'text' 'tag' 'lbl_lowpass' 'string' 'Cutoff (Hz):' 'horizontalalignment' 'left' 'enable' onoff_lp});
add({'style' 'edit' 'tag' 'lowpass' 'string' iff(isempty(choices.lowpass),'',num2str(choices.lowpass)) 'enable' onoff_lp});

% Segmentation
add({'style' 'text' 'string' 'Segmentation (epoching):' 'horizontalalignment' 'left'});
add({'style' 'checkbox' 'tag' 'apply_epoch' 'value' choices.apply_epoch 'string' 'Enable' 'callback' cb_seg});
add({'style' 'text' 'string' ''});
add({'style' 'text' 'tag' 'lbl_epoch' 'string' 'Epoch window [ms] (start end):' 'horizontalalignment' 'left' 'enable' onoff_seg});
add({'style' 'edit' 'tag' 'epoch_window' 'string' sprintf('%d %d',choices.epoch_window) 'enable' onoff_seg});

% CAR
add({'style' 'text' 'string' 'Common Average Reference (CAR):' 'horizontalalignment' 'left'});
add({'style' 'checkbox' 'tag' 'apply_acar' 'value' choices.apply_acar 'string' 'Enable' 'callback' cb_acar});
add({'style' 'text' 'string' ''});
add({'style' 'text' 'tag' 'lbl_acar_fraction' 'string' 'Fraction of channels (0–1):' 'horizontalalignment' 'left' 'enable' onoff_car});
add({'style' 'edit' 'tag' 'acar_fraction' 'string' num2str(choices.acar_fraction) 'enable' onoff_car});
add({'style' 'text' 'string' ''});
add({'style' 'text' 'tag' 'lbl_acar_timewin' 'string' 'CAR window [ms] (start end):' 'horizontalalignment' 'left' 'enable' onoff_car});
add({'style' 'edit' 'tag' 'acar_timewin' 'string' sprintf('%d %d',choices.acar_timewin) 'enable' onoff_car});

% Baseline
add({'style' 'text' 'string' 'Remove baseline:' 'horizontalalignment' 'left'});
add({'style' 'checkbox' 'tag' 'apply_baseline' 'value' choices.apply_baseline 'string' 'Enable' 'callback' cb_bl});
add({'style' 'text' 'string' ''});
add({'style' 'text' 'tag' 'lbl_bl_method' 'string' 'Method:' 'horizontalalignment' 'left' 'enable' onoff_bl});
add({'style' 'popupmenu' 'tag' 'bl_method' 'string' {'Median (default)' 'Mean' 'Trimmed mean' '1/f'} 'value' choices.baseline_method 'enable' onoff_bl});
add({'style' 'text' 'string' ''});
add({'style' 'text' 'tag' 'lbl_bl_period' 'string' 'Baseline period [ms] (start end):' 'horizontalalignment' 'left' 'enable' onoff_bl});
add({'style' 'edit' 'tag' 'bl_period' 'string' sprintf('%d %d',choices.baseline_period) 'enable' onoff_bl});
add({'style' 'text' 'string' ''});
add({'style' 'text' 'tag' 'lbl_bl_mode' 'string' 'Mode:' 'horizontalalignment' 'left' 'enable' onoff_bl});
add({'style' 'popupmenu' 'tag' 'bl_mode' 'string' {'Subtract (default)' 'Divide'} 'value' choices.baseline_mode 'enable' onoff_bl});

% launch (narrow)
[res, ~, ~, out] = inputgui('geometry', uigeom, 'uilist', uilist, ...
                            'title', 'iEEGLAB: Preprocessing Options', ...
                            'minwidth', 440);
% Clean abort on Cancel or window close
if isempty(res) || isempty(out)
    fprintf('iEEGLAB preprocess dialog canceled\n');
    wasCanceled = true;  
    return
end

% parse outputs
if isfield(out,'evsel_json') && ~isempty(out.evsel_json)
    try, choices.event_filters = jsondecode(out.evsel_json); catch, choices.event_filters = struct(); end
end

if isfield(out,'min_trials')
    mt = str2double(out.min_trials); if isfinite(mt) && mt>=0, choices.min_trials = mt; end
end
choices.remove_rare_cond = isfield(out,'remove_rare_cond') && logical(out.remove_rare_cond);
choices.remove_no_coords = isfield(out,'rm_no_coords') && logical(out.rm_no_coords);

choices.apply_ds = isfield(out,'apply_ds') && logical(out.apply_ds);
if isfield(out,'ds_rate') && ~isempty(out.ds_rate)
    dr = str2double(out.ds_rate); if isfinite(dr) && dr>0, choices.downsample = dr; end
end

choices.apply_highpass = isfield(out,'apply_highpass') && logical(out.apply_highpass);
if isfield(out,'highpass') && ~isempty(out.highpass)
    hpv = str2double(out.highpass); if isfinite(hpv) && hpv>0, choices.highpass = hpv; end
end
if isfield(out,'highpass_type') && ~isempty(out.highpass_type)
    choices.highpass_type = out.highpass_type;
end
hp_labels = {'noncausal-zero-phase','causal-minimum-phase'};
choices.highpass_type_label = hp_labels{min(max(choices.highpass_type,1),2)};

choices.apply_notch = isfield(out,'apply_notch') && logical(out.apply_notch);
if isfield(out,'notch') && ~isempty(out.notch)
    nv = str2double(out.notch); if isfinite(nv) && nv>0, choices.notch = nv; end
end

choices.apply_lowpass = isfield(out,'apply_lowpass') && logical(out.apply_lowpass);
if isfield(out,'lowpass') && ~isempty(out.lowpass)
    lp = str2double(out.lowpass);
    if isfinite(lp) && lp>0, choices.lowpass = lp; else, choices.lowpass = []; end
else
    if ~choices.apply_lowpass, choices.lowpass = []; end
end

choices.apply_epoch = isfield(out,'apply_epoch') && logical(out.apply_epoch);
if isfield(out,'epoch_window') && ~isempty(out.epoch_window)
    tw = sscanf(out.epoch_window,'%f'); if numel(tw)>=2, choices.epoch_window = tw(1:2).'; end
end

choices.apply_acar = isfield(out,'apply_acar') && logical(out.apply_acar);
if isfield(out,'acar_fraction') && ~isempty(out.acar_fraction)
    choices.acar_fraction = max(0, min(1, str2double(out.acar_fraction)));
end
if isfield(out,'acar_timewin') && ~isempty(out.acar_timewin)
    tw = sscanf(out.acar_timewin,'%f'); if numel(tw)>=2, choices.acar_timewin = tw(1:2).'; end
end

choices.apply_baseline = isfield(out,'apply_baseline') && logical(out.apply_baseline);
bl_labels = {'median','mean','trimmed mean','1/f'};
if isfield(out,'bl_method') && ~isempty(out.bl_method), bl_idx = out.bl_method; else, bl_idx = choices.baseline_method; end
choices.baseline_method = bl_labels{min(max(bl_idx,1),4)};
if isfield(out,'bl_period') && ~isempty(out.bl_period)
    bp = sscanf(out.bl_period,'%f'); if numel(bp)>=2, choices.baseline_period = bp(1:2).'; end
end
bm_labels = {'subtract','divide'};
if isfield(out,'bl_mode') && ~isempty(out.bl_mode), bm_idx = out.bl_mode; else, bm_idx = choices.baseline_mode; end
choices.baseline_mode = bm_labels{min(max(bm_idx,1),2)};

%% merge into EEG.ieeglab.opt

if ~isfield(EEG,'ieeglab'), EEG.ieeglab = struct(); end
if ~isfield(EEG.ieeglab,'opt') || ~isstruct(EEG.ieeglab.opt), EEG.ieeglab.opt = struct(); end
EEG.ieeglab.opt = merge_structs(EEG.ieeglab.opt, choices);

end


%% sub-GUI (unchanged)
function local_open_event_selector(srcBtn, evT, showCols)
lists = cell(1,numel(showCols)); tags = cell(1,numel(showCols));
for k = 1:numel(showCols)
    col = showCols{k}; vals = evT.(col);
    if iscell(vals), vals = string(vals); end
    if iscategorical(vals), vals = string(vals); end
    if isnumeric(vals), vals = string(vals); end
    u = unique(vals(~ismissing(vals)));
    lists{k} = [{'All'}; cellstr(u(:))];
    tags{k}  = ['fld_' regexprep(col,'[^\w]','_')];
end
nRows = 1+numel(showCols); uigeom = cell(nRows,1); gvert = ones(1,nRows);
uilist = cell(0,1);
uilist{end+1} = {'style' 'text' 'string' 'Select events to KEEP (leave at "All")' 'fontweight' 'bold'};
uigeom{1} = [1]; gvert(1)=1; ii=2;
for k=1:numel(showCols)
    uigeom{ii}=[1 1]; gvert(ii)=1.4;
    uilist{end+1}={'style' 'text' 'string' showCols{k} 'horizontalalignment' 'left'};
    uilist{end+1}={'style' 'listbox' 'tag' tags{k} 'string' lists{k} 'value' 1 'min' 0 'max' 2};
    ii=ii+1;
end
[res,~,~,out]=inputgui('geometry',uigeom,'geomvert',gvert,'uilist',uilist,'title','Select events of interest','minwidth',440);
if isempty(res)||isempty(out), return; end
filters=struct();
for k=1:numel(showCols)
    tg=tags{k}; selIdx=1;
    if isfield(out,tg)&&~isempty(out.(tg)), selIdx=out.(tg); end
    if any(selIdx==1), filters.(showCols{k})={}; else, filters.(showCols{k})=lists{k}(selIdx); end
end
hFig=ancestor(srcBtn,'figure'); hStore=findobj(hFig,'tag','evsel_json');
if ~isempty(hStore)&&ishghandle(hStore)
    try, set(hStore,'string',jsonencode(filters));
    catch
        keys=fieldnames(filters); fr=strings(0);
        for i=1:numel(keys)
            v=string(filters.(keys{i})); if isempty(v), v="<ALL>"; end
            fr(end+1)=keys{i}+"="+strjoin(v,"|"); %#ok<AGROW>
        end
        set(hStore,'string',char(strjoin(fr,";")));
    end
end
end

%% helpers

function out = iff(cond,a,b)
    if cond, out=a; else, out=b; end
end

function out = tern(x)
    if x, out='on'; else, out='off'; end
end

function s = merge_structs(a,b)
    s=a; f=fieldnames(b);
    for i=1:numel(f), s.(f{i})=b.(f{i}); end
end

function cb = local_cb_toggle(tag_list)
    % Returns a function handle callback that toggles enable on controls with given tags
    cb = @(h,~) local_set_enable(ancestor(h,'figure'), tag_list, get(h,'value')~=0);
end

function local_set_enable(hFig, tag_list, is_on)
    onoff = ternary(is_on, 'on', 'off');
    for k = 1:numel(tag_list)
        h = findobj(hFig,'tag',tag_list{k});
        if ~isempty(h) && ishghandle(h)
            set(h,'enable',onoff);
        end
    end
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
