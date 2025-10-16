function signal_clean = reject_bad_trials_ccep(signal, stim_labels, chan_idx, k_thresh)
% Reject bad trials in CCEP data relative to stim current level
%
% signal      - chan x time x trials
% stim_labels - cell array {trials} e.g., {'4.0 mA', '4.0 mA', '6.0 mA', ...}
% chan_idx    - channel index or vector of channels to evaluate amplitude
% k_thresh    - threshold multiplier for MAD (default 6)
%
% Returns:
% signal_clean - signal with bad trials removed

if nargin < 4 || isempty(k_thresh), k_thresh = 6; end

[nCh, nTime, nTrials] = size(signal);
if length(stim_labels) ~= nTrials
    error('Length of stim_labels (%d) must match number of trials (%d)', ...
        length(stim_labels), nTrials);
end

kept = true(1, nTrials);

stim_levels = unique(stim_labels);
for s = 1:numel(stim_levels)
    idx = strcmp(stim_labels, stim_levels{s});
    sig_s = squeeze(signal(chan_idx, :, idx));

    % Peak-to-peak per trial
    if ndims(sig_s) == 3
        p2p = squeeze(max(sig_s,[],2) - min(sig_s,[],2)); % chan x trials
        p2p = max(p2p, [], 1); % max across channels
    else
        p2p = max(sig_s,[],2) - min(sig_s,[],2); % trials
    end

    % Robust stats
    med_val = median(p2p);
    mad_val = mad(p2p,1);

    % Avoid divide-by-zero if MAD is 0
    if mad_val == 0
        bad_idx = false(size(p2p));
    else
        bad_idx = p2p > (med_val + k_thresh * mad_val);
    end

    kept(idx) = ~bad_idx;
end

fprintf('Rejected %d/%d trials (%.1f%%)\n', sum(~kept), nTrials, 100*mean(~kept));

% Remove bad trials
signal_clean = signal(:,:,kept);

end
