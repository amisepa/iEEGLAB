function [opt, wasCanceled] = ieeglab_gui_load1(filepath)
% ieeglab_gui_load1
%
% Open a small GUI to collect iEEG loading parameters (TSV paths only).
% Cleanly aborts on Cancel or window close.
%
% Inputs (optional)
%   filepath     char/string: default folder to open file choosers in
%
% Outputs
%   opt          struct with fields (empty [] when canceled):
%                  .elec_tsv
%                  .events_tsv
%   wasCanceled  logical, true if user canceled

wasCanceled = false;      % default

% Resolve start_dir from input or fallbacks
if nargin >= 1 && ~isempty(filepath) && isfolder(filepath)
    start_dir = char(filepath);
else
    plugin_path = fileparts(which('eegplugin_ieeglab'));
    if isempty(plugin_path)
        start_dir = pwd;
    else
        start_dir = fullfile(plugin_path, 'tutorial');
        if ~isfolder(start_dir)
            start_dir = plugin_path;
        end
    end
end

% --- Callbacks (function handles; fewer quoting issues) ---
    function do_browse_elec(src, ~)
        fprintf('\nSelect the electrodes.tsv file...\n');
        [f,p] = uigetfile({'*.tsv','Electrode TSV (*.tsv)'}, ...
                          'Select electrode locations file', start_dir);
        if isequal(f,0)
            fprintf('Canceled selecting electrode locations\n');
            return;
        end
        set(findobj(ancestor(src,'figure'),'tag','elec_tsv'), 'string', fullfile(p,f));
        hFig = ancestor(src,'figure'); try, figure(hFig); uistack(hFig,'top'); drawnow; end
    end

    function do_browse_events(src, ~)
        fprintf('\nSelect the events.tsv file...\n');
        [f,p] = uigetfile({'*.tsv','Events TSV (*.tsv)'}, ...
                          'Select events file', start_dir);
        if isequal(f,0)
            fprintf('Canceled selecting events\n');
            return;
        end
        set(findobj(ancestor(src,'figure'),'tag','events_tsv'), 'string', fullfile(p,f));
        hFig = ancestor(src,'figure'); try, figure(hFig); uistack(hFig,'top'); drawnow; end
    end
% ----------------------------------------------------------

% GUI layout
uilist = {
    {'style' 'text'  'string' 'Electrode locations file (.tsv):' 'fontweight' 'bold'}
    {'style' 'edit'  'string' '' 'tag' 'elec_tsv' 'enable' 'inactive' 'horizontalalignment' 'left'}
    {'style' 'pushbutton' 'string' 'Browse…' 'callback' @do_browse_elec}
    {}
    {'style' 'text'  'string' 'Events file (.tsv):' 'fontweight' 'bold'}
    {'style' 'edit'  'string' '' 'tag' 'events_tsv' 'enable' 'inactive' 'horizontalalignment' 'left'}
    {'style' 'pushbutton' 'string' 'Browse…' 'callback' @do_browse_events}
};

uigeom = {
    [.45 .4 .22]
    1
    [.45 .4 .22]
};

% Launch GUI
[res, ~, ~, outstruct] = inputgui(uigeom, uilist, 'pophelp(''ieeglab_load'')', ...
                                  'iEEGLAB plugin  Loading (TSV only)');

% Clean abort on Cancel or window close
if isempty(res) || isempty(outstruct)
    fprintf('iEEGLAB load dialog canceled\n');
    wasCanceled = true;
    opt = [];
    return;
end

% Map results from uilist controls that produce values (order matches uilist edits)
elec_tsv_value   = char(res{1});
events_tsv_value = char(res{2});

% Basic validation (warnings only)
if ~isempty(elec_tsv_value) && ~endsWith(lower(elec_tsv_value), '.tsv')
    warning('Electrode file does not have .tsv extension: %s', elec_tsv_value);
end
if ~isempty(events_tsv_value) && ~endsWith(lower(events_tsv_value), '.tsv')
    warning('Events file does not have .tsv extension: %s', events_tsv_value);
end

opt = struct();
opt.elec_tsv   = elec_tsv_value;
opt.events_tsv = events_tsv_value;

fprintf('\nSelected parameters\n');
disp(opt);
end
