function [signal, kept] = read_mefd_data(mefd_base_dir, channel_names, start_smp, stop_smp)
% Read MEF 3.0 data trials for multiple channels
%
% INPUTS:
%   mefd_base_dir  - directory containing channel folders (e.g., 'C3.timd', 'C4.timd', etc.)
%   channel_names  - cell array of channel names (matching folder names without extension)
%   start_smp      - vector of trial start samples
%   stop_smp       - vector of trial stop samples
%
% OUTPUTS:
%   signal - [nChan x nTrials x nTime] array of trial data
%   kept   - logical vector of trials successfully read

    
    % Create MEF3 object (empty password if not encrypted)
    % MEF3 = MEF_3p0.MultiscaleElectrophysiologyFile_3p0();
    MEF3 = MultiscaleElectrophysiologyFile_3p0();

    nChan   = numel(channel_names);
    nTrials = numel(start_smp);
    nTime   = stop_smp(1) - start_smp(1) + 1;

    signal = NaN(nChan, nTrials, nTime);
    kept   = false(nTrials, 1);

    for t = 1:nTrials
        for c = 1:nChan
            try
                chan_dir = fullfile(mefd_base_dir, [channel_names{c} '.timd']);
                % Call MEF 3.0 reader
                data = read_mef_data(MEF3, chan_dir, '', 'samples', start_smp(t), stop_smp(t));
                signal(c,t,:) = data;
            catch ME
                warning('Failed to read chan %s trial %d: %s', channel_names{c}, t, ME.message);
                signal(c,t,:) = NaN;
            end
        end
        kept(t) = true;
    end
end
