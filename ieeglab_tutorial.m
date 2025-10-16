%% Demo1 - CCEP
% Tutorial dataset is subject 01 from BIDS dataset available here: https://openneuro.org/datasets/ds004977/versions/1.2.0
% downsampled to 100 Hz a few channels were preserved.
% 
% Copyright (C) - iEEGLAB Team, 2025

clear; close all; clc

plugin_path = find(fileparts(which()));

eeglab; close


%% Demo2 (from current eeglab tutorial)

% https://dandiarchive.org/dandiset/000576/
% https://nemar.org/dataexplorer/detail?dataset_id=ds003708&processed=0