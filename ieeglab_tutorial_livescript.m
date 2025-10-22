%[text] # Step 1: launch EEGLAB
clear; close all; clc
eeglab; close %[output:7188b9da]

plugin_path = fileparts(which('eegplugin_ieeglab'));
cd(plugin_path)
%%
%[text] ## Load eCoG .vhdr data

filepath = '/Users/cedriccannard/Documents/ecog/sub-02/ses-01/ieeg/';
filename = 'sub-02_ses-01_task-visual_run-01_ieeg.vhdr';
EEG = pop_loadbv(filepath,filename);
EEG.filepath = filepath;
EEG.filename = 'eCoG dataset';
EEG.setname  = EEG.filename;
EEG.srate    = round(EEG.srate);
EEG = eeg_checkset(EEG);
%%
%[text] ## Load long sEEG .mefd dataset
filepath = '/Users/cedriccannard/Documents/dataset3';
filename = 'sub-02_ses-ieeg01_task-ccep_run-01_ieeg.mefd';
% EEG = pop_MEF3(fullfile(filepath, filename));


% [metadata, signal] = ieeglab_load_mefd(fileName);
metadata     = readMef3(fullfile(filepath, filename));
fs           = double(metadata.time_series_metadata.section_2.sampling_frequency);
startSample  = 0; % 369
endSample    = 3000 * fs;
% parpool;
tic
[metadata, signal] = ieeglab_load_mefd(fullfile(filepath, filename), [], [], ...
    'samples', [startSample endSample]);
toc %[output:4478379f]

% Convert to EEGLAB format
EEG = eeg_emptyset;
EEG.filepath = filepath;
EEG.data = signal;
EEG.nbchan = size(EEG.data,1);
EEG.pnts = size(EEG.data,2);
EEG.srate = fs;
EEG = eeg_checkset(EEG); %[output:1958d0e9]

% % Channel labels from metadata
% for iChan = 1:EEG.nbchan
%     EEG.chanlocs(iChan).labels = metadata.time_series_channels(iChan).name;
% end
% EEG = eeg_checkset(EEG);

%%
%[text] ## sEEG .mefd data1 (Dora's reduced dataset)
filepath = '/Users/cedriccannard/Documents/dataset1';
filename = 'sub-01_ses-ieeg01_task-ccep_run-01_ieeg.mefd';
EEG = pop_MEF3(fullfile(filepath, filename));
EEG.filepath = filepath;

%%
%[text] ## Load 3D Cartesian (XYZ) electrode locations & events
EEG = ieeglab_load(EEG); %[output:79823838]

%%
%[text] ## Visualize
% Plot the raw time series
pop_eegplot(EEG,1,1,1);


% Electrodes in 3D glass brain 
addpath(genpath('/Users/cedriccannard/Documents/MATLAB/vistasoft'))
% addpath('/Users/cedriccannard/Documents/MATLAB/fieldtrip')
ieeglab_vis_elec(EEG);

%%
%[text] ## Preprocess

EEG = ieeglab_preprocess(EEG);

%%
%[text] ## Some classic EEGLAB plots


for iChan = 1:25:EEG.nbchan
    figure
    pop_erpimage(EEG,1, iChan,[],EEG.chanlocs(iChan).labels,10,1,{},[],'' ,'yerplabel','\muV','erp','on','cbar','on');
    pause(1)
    close(gcf)
end

figure; 
pop_timtopo(EEG, [EEG.times(1) EEG.times(end)], [], ...
    'CCEP - all electrodes','verbose','off');

% --> EDIT TO SHOW THE 3D GLASS BRAIN INSTEAD
% figure; pop_timtopo(EEG, [-500  990], NaN, 'ERP data and scalp maps');

% figure; pop_plottopo(EEG, 1:EEG.nbchan , '', 0, 'ydir',1);
figure; plottopo( trimmean(EEG.data,20,3), 'frames', EEG.pnts, 'limits', [-500 990 0 0], 'chans', 1:EEG.nbchan, 'ydir', 1);
% figure; plottopo( EEG.data, 'frames', EEG.pnts, 'limits', [-500 990 0 0], 'chans', 1:EEG.nbchan, 'ydir', 1);

% PSD
% [pwr, pwr_osc, psd, psd_osc, f] = compute_pwr(EEG.data, EEG.srate,'PlotPSD', true, 'UseParallel', true);
figure; pop_spectopo(EEG, 1, [], 'EEG' , 'freq', [], 'freqrange',[0.1 50],'electrodes','off');

% Heatmap of the average across trials (all channels x time)
avg_pair = squeeze(trimmean(EEG.data, 10, 3));             % [chan x time]
plot_ccep(avg_pair, EEG.times, {EEG.event.type}, 'all', [], 0.20);

% single-channel overlay of all trials + mean
channel = 10;   
trials = contains({EEG.event.type}, EEG.chanlocs(channel).labels);
sum(trials)
plot_ccep(EEG.data(:,:,trials), EEG.times, {EEG.chanlocs.labels}, 'single', channel, []);



%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
%[output:7188b9da]
%   data: {"dataType":"text","outputData":{"text":"Disable MATLAB Copilot in settings to hide \"Explain Error\" buttons\nSome menu items hidden. Use Preference menu to show them all.\neeglab: options file is ~\/eeg_options.m\nRetrieving plugin versions from server...\nRetrieving download statistics...\nEEGLAB: adding \"AMICA\" v1.7 (see >> help eegplugin_amica)\nEEGLAB: adding \"Biosig\" v3.8.4 to the path\nEEGLAB: adding \"EEG-BIDS\" v10.3 (see >> help eegplugin_eegbids)\nEEGLAB: adding \"ICLabel\" v1.7 (see >> help eegplugin_iclabel)\nEEGLAB: adding \"MEF3v\" v1.2.2 (see >> help eegplugin_MEF3)\nEEGLAB: adding \"bva-io\" v1.74 (see >> help eegplugin_bva_io)\nEEGLAB: adding \"clean_rawdata\" v2.11 (see >> help eegplugin_clean_rawdata)\nEEGLAB: adding \"dipfit\" v5.6 (see >> help eegplugin_dipfit)\nEEGLAB: adding \"firfilt\" v2.8 (see >> help eegplugin_firfilt)\nEEGLAB: adding \"iEEGLAB\" v1.0 (see >> help eegplugin_ieeglab)\nEEGLAB: adding \"neuroscanio\" v1.8 (see >> help eegplugin_neuroscanio)\nEEGLAB: adding \"nwbio\" v1.2 (see >> help eegplugin_nwbio)\nEEGLAB: adding \"roiconnect\" v1.1 (see >> help eegplugin_roiconnect)\nDisable MATLAB Copilot in settings to hide \"Explain Error\" buttons\nYou are using the latest version of EEGLAB.\n","truncated":false}}
%---
%[output:4478379f]
%   data: {"dataType":"text","outputData":{"text":"Elapsed time is 548.899495 seconds.\n","truncated":false}}
%---
%[output:1958d0e9]
%   data: {"dataType":"text","outputData":{"text":"eeg_checkset warning: 3rd dimension size of data (1) does not match the number of epochs (0), corrected\neeg_checkset note: upper time limit (xmax) adjusted so (xmax-xmin)*srate+1 = number of frames\n","truncated":false}}
%---
%[output:79823838]
%   data: {"dataType":"text","outputData":{"text":"iEEGLAB load dialog canceled\n","truncated":false}}
%---
