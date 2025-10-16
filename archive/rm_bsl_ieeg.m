function signal = rm_bsl_ieeg(signal, bsl_idx, method, fs)
% rm_bsl_ieeg - Remove baseline from iEEG/EEG data epochs
%
% Usage:
%   signal = rm_bsl_ieeg(signal, bsl_idx)                       % default = median
%   signal = rm_bsl_ieeg(signal, bsl_idx, method)
%   signal = rm_bsl_ieeg(signal, bsl_idx, '1/f', fs)
%
% Inputs:
%   signal   - Data matrix [nChannels x nTimePoints x nEpochs]
%   bsl_idx  - Indices for baseline calculation (vector of sample indices)
%   method   - 'median' (default), 'mean', 'trimmed_mean', or '1/f'
%              'mean'         : subtract mean of baseline window
%              'median'       : subtract median of baseline window
%              'trimmed_mean' : subtract 20%% trimmed mean of baseline window
%              '1/f'          : fit broadband 1/f slope+offset (additive model)
%                               to baseline spectrum (excluding oscillatory bands)
%                               and subtract from entire trial
%   fs       - Sampling rate in Hz (required only for '1/f' method)
%
% Outputs:
%   signal   - Baseline-corrected data
%
% Cedric Cannard © iEEGLAB Plugin, 2025
%

%% --- Defaults and validation ---
if nargin < 3 || isempty(method)
    method = 'median';
else
    method = lower(method);
end

valid_methods = {'mean','median','trimmed_mean','1/f'};
if ~ismember(method, valid_methods)
    error('Invalid method. Choose: ''mean'', ''median'', ''trimmed_mean'', or ''1/f''.');
end

if ndims(signal) ~= 3
    error('signal must be [nChannels x nTimePoints x nEpochs]');
end
[nCh, nTim, nEp] = size(signal);

if any(bsl_idx < 1) || any(bsl_idx > nTim)
    error('Baseline indices must be within 1–%d', nTim);
end
if isempty(bsl_idx)
    warning('Empty baseline indices: no correction applied.');
    return;
end

if strcmp(method, '1/f')
    if nargin < 4 || isempty(fs) || ~isscalar(fs) || fs <= 0
        error('Sampling rate fs (Hz) is required for ''1/f'' baseline removal.');
    end
end

%% --- Status message ---
fprintf('Removing baseline using "%s" method...\n', method);

if strcmp(method, '1/f')
    fprintf(['[INFO] Using additive-model 1/f correction.\n' ...
        '       Fits and subtracts broadband slope/offset from baseline.\n' ...
        '       For iEEG, this is mainly relevant if broadband activity\n' ...
        '       varies across conditions, electrodes, or brain states.\n']);
end

%% --- Baseline correction ---
switch method
    case 'mean'
        bsl_vals = mean(signal(:,bsl_idx,:), 2); % mean over baseline timepoints
        signal = signal - repmat(bsl_vals, [1 nTim 1]);

    case 'median'
        bsl_vals = median(signal(:,bsl_idx,:), 2);
        signal = signal - repmat(bsl_vals, [1 nTim 1]);

    case 'trimmed_mean'
        trim_pct = 20; % 20% trimming
        trim_amt = round((trim_pct/100) * numel(bsl_idx));
        bsl_vals = nan(nCh, 1, nEp);
        for ch = 1:nCh
            for ep = 1:nEp
                vals = squeeze(signal(ch,bsl_idx,ep));
                vals = sort(vals);
                if numel(vals) > 2*trim_amt
                    vals = vals(trim_amt+1:end-trim_amt);
                end
                bsl_vals(ch,1,ep) = mean(vals);
            end
        end
        signal = signal - repmat(bsl_vals, [1 nTim 1]);

    case '1/f'
        % Frequency vector
        f = (0:nTim-1)' * (fs/nTim);
        exclude_band = [4 20]; % Hz to exclude from fit
        fit_mask = f > 0 & (f < exclude_band(1) | f > exclude_band(2)) & f <= fs/2;

        for ch = 1:nCh
            for ep = 1:nEp
                % Ensure column vector for FFT
                x = squeeze(signal(ch,:,ep)).'; % nTim x 1 column
                xb = x(bsl_idx);

                % FFT of baseline
                Pb = abs(fft(xb, nTim)).^2;

                if sum(fit_mask) < 5
                    warning('Too few points for 1/f fit (ch=%d, ep=%d)', ch, ep);
                    continue;
                end

                % Fit slope/intercept in log-log space
                f_log = log10(f(fit_mask));
                P_log = log10(Pb(fit_mask));
                p_fit = polyfit(f_log, P_log, 1);
                slope = p_fit(1);
                intercept = p_fit(2);

                % Predict broadband spectrum for full trial
                P1f_pred = 10.^(intercept + slope * log10(f + eps));

                % FFT of full trial
                X_full = fft(x, nTim);
                mag = abs(X_full);
                phase = angle(X_full);

                % Subtract broadband magnitude
                mag_corrected = mag - sqrt(P1f_pred(:));
                mag_corrected(mag_corrected < 0) = 0;

                % Reconstruct time series
                X_corrected = mag_corrected .* exp(1i * phase);
                x_corrected = real(ifft(X_corrected, nTim));

                % Assign back as 1 x nTim row vector
                signal(ch,:,ep) = reshape(x_corrected, 1, nTim);
            end
        end
end

