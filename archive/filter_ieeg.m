function signal_filt = filter_ieeg(signal, fs, varargin)
% filter_ieeg - Apply EEGLAB's pop_eegfiltnew FIR filter to epoched data.
%
% Usage:
%   signal_filt = filter_ieeg(signal, fs, 'locutoff', 0.5, 'hicutoff', [], ...
%                              'revfilt', 0, 'usefftfilt', 0, 'plotfreqz', 0)
%
% Inputs:
%   signal     - [nChan x nTime x nTrials] data matrix
%   fs         - Sampling rate (Hz)
%
% Optional name/value pairs (directly passed to pop_eegfiltnew):
%   'locutoff'   - Lower passband edge in Hz ([] or 0 → lowpass)
%   'hicutoff'   - Upper passband edge in Hz ([] or 0 → highpass)
%   'filtorder'  - Filter order (even number, default auto)
%   'revfilt'    - [0|1] invert filter (bandpass→notch)
%   'plotfreqz'  - [0|1] plot freq/phase response
%   'minphase'   - [true|false] min-phase causal filter
%   'usefftfilt' - [0|1] use FFT-based convolution
%
% Output:
%   signal_filt - Filtered data, same size as input
%
% Notes:
%   - Requires EEGLAB on MATLAB path.
%   - This wraps pop_eegfiltnew() for 3D epoched arrays.
%
% Cedric Cannard © iEEGLAB plugin, 2025

% --- Check EEGLAB availability ---
if ~exist('pop_eegfiltnew', 'file')
    error('EEGLAB function pop_eegfiltnew() not found. Add EEGLAB to your MATLAB path.');
end

% --- Get input size ---
[nChan, nTime, nTrials] = size(signal);

% --- Create temporary EEGLAB EEG struct ---
EEG = eeg_emptyset();
EEG.data    = reshape(signal, nChan, nTime * nTrials);
EEG.srate   = fs;
EEG.nbchan  = nChan;
EEG.trials  = nTrials;
EEG.pnts    = nTime;
EEG.xmin    = 0;
EEG.xmax    = (nTime-1)/fs;

% --- Filter using EEGLAB's pop_eegfiltnew ---
EEG = pop_eegfiltnew(EEG, varargin{:});

% --- Reshape back to [chan × time × trials] ---
signal_filt = reshape(EEG.data, nChan, nTime, nTrials);
end
