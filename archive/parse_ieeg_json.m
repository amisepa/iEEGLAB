function info = parse_ieeg_json(json_path)
% Parse BIDS iEEG JSON sidecar for key recording parameters.
%
% INPUT:
%   json_path - full path to the JSON file (string or char)
%
% OUTPUT:
%   info - struct with fields:
%       .SamplingFrequency  (Hz)
%       .PowerLineFrequency (Hz)
%       .RecordingDuration  (sec)
%       .HighpassCutoff     (Hz)
%       .LowpassCutoff      (Hz)
%       .StimParameters     (string)

    if ~exist(json_path, 'file')
        error('JSON file not found: %s', json_path);
    end

    % Read JSON
    raw = fileread(json_path);
    data = jsondecode(raw);

    info = struct();

    % --- Required ---
    if isfield(data, 'SamplingFrequency') && ~isempty(data.SamplingFrequency)
        info.SamplingFrequency = data.SamplingFrequency;
        fprintf('Sampling frequency: %g Hz\n', info.SamplingFrequency);
    else
        error('Required field "SamplingFrequency" not found in %s', json_path);
    end

    % --- Optional ---
    if isfield(data, 'PowerLineFrequency') && ~isempty(data.PowerLineFrequency)
        info.PowerLineFrequency = data.PowerLineFrequency;
        fprintf('Power line frequency: %g Hz\n', info.PowerLineFrequency);
    else
        warning('PowerLineFrequency not found in JSON.');
        info.PowerLineFrequency = [];
    end

    if isfield(data, 'RecordingDuration') && ~isempty(data.RecordingDuration)
        info.RecordingDuration = data.RecordingDuration;
        fprintf('Recording duration: %g seconds\n', info.RecordingDuration);
    else
        warning('RecordingDuration not found in JSON.');
        info.RecordingDuration = [];
    end

    % --- Hardware filters ---
    if isfield(data, 'HardwareFilters') && isstruct(data.HardwareFilters)
        if isfield(data.HardwareFilters, 'HighpassFilter') && ...
           isfield(data.HardwareFilters.HighpassFilter, 'CutoffFrequency')
            info.HighpassCutoff = data.HardwareFilters.HighpassFilter.CutoffFrequency;
            fprintf('Highpass filter cutoff: %g Hz\n', info.HighpassCutoff);
        else
            warning('Highpass filter cutoff not found.');
            info.HighpassCutoff = [];
        end

        if isfield(data.HardwareFilters, 'LowpassFilter') && ...
           isfield(data.HardwareFilters.LowpassFilter, 'CutoffFrequency')
            info.LowpassCutoff = data.HardwareFilters.LowpassFilter.CutoffFrequency;
            fprintf('Lowpass filter cutoff: %g Hz\n', info.LowpassCutoff);
        else
            warning('Lowpass filter cutoff not found.');
            info.LowpassCutoff = [];
        end
    else
        warning('HardwareFilters struct not found.');
        info.HighpassCutoff = [];
        info.LowpassCutoff = [];
    end

    % --- Stimulation parameters ---
    if isfield(data, 'ElectricalStimulationParameters') && ~isempty(data.ElectricalStimulationParameters)
        info.StimParameters = data.ElectricalStimulationParameters;
        fprintf('Stimulation parameters: %s\n', info.StimParameters);
    else
        warning('ElectricalStimulationParameters not found.');
        info.StimParameters = '';
    end
end
