function EEG = ieeglab_rm_baseline(EEG)
% ieeglab_rm_baseline  Apply baseline correction to EEG.data using GUI choices.
%
% Uses:
%   EEG.data                : [channels x samples x trials]
%   EEG.srate               : sampling rate (Hz) for '1/f' method
%   EEG.times               : time (ms), length must match samples
%   EEG.ieeglab.opt fields  :
%       baseline_method     : 'median' | 'mean' | 'trimmed mean' | '1/f'
%       baseline_period     : [start_ms end_ms]
%       baseline_mode       : 'subtract' | 'divide'
%
% Output:
%   EEG.data (baseline-corrected)

% ---------- sanity ----------
assert(isfield(EEG,'data') && ~isempty(EEG.data), 'EEG.data is empty.');
X = EEG.data;
if ndims(X)==2, X = reshape(X, size(X,1), size(X,2), 1); end
[C,T,N] = size(X);

assert(isfield(EEG,'times') && numel(EEG.times)==T, ...
    'EEG.times must exist and match the time dimension of EEG.data.');

% ---------- options from GUI ----------
opt = struct('baseline_method','median', 'baseline_period',[-500 -50], 'baseline_mode','subtract');
if isfield(EEG,'ieeglab') && isfield(EEG.ieeglab,'opt') && isstruct(EEG.ieeglab.opt)
    f = EEG.ieeglab.opt;
    if isfield(f,'baseline_method') && ~isempty(f.baseline_method)
        opt.baseline_method = f.baseline_method;
    end
    if isfield(f,'baseline_period') && numel(f.baseline_period)>=2
        opt.baseline_period = double(f.baseline_period(1:2));
    end
    if isfield(f,'baseline_mode') && ~isempty(f.baseline_mode)
        opt.baseline_mode = f.baseline_mode;
    end
end

% Normalize names
method = lower(strtrim(string(opt.baseline_method)));
method = strrep(method,'_',' ');                 % allow 'trimmed_mean'
mode   = lower(strtrim(string(opt.baseline_mode)));
valid_methods = ["median","mean","trimmed mean","1/f"];
assert(any(method==valid_methods), 'Invalid baseline_method: %s', method);
assert(any(mode==["subtract","divide"]), 'Invalid baseline_mode: %s', mode);

% '1/f' needs fs
if method=="1/f"
    assert(isfield(EEG,'srate') && isscalar(EEG.srate) && EEG.srate>0, ...
        'EEG.srate required for 1/f baseline.');
    if mode=="divide"
        warning('Mode "divide" not supported for 1/f; using "subtract" instead.');
        mode = "subtract";
    end
end

% ---------- baseline indices from EEG.times (ms) ----------
t_ms = double(EEG.times(:));
bsl_ms = sort(double(opt.baseline_period(:))).';
[~,i0] = min(abs(t_ms - bsl_ms(1)));
[~,i1] = min(abs(t_ms - bsl_ms(2)));
i0 = max(1, min(T, i0));
i1 = max(1, min(T, i1));
if i1 < i0, [i0,i1] = deal(i1,i0); end
bsl_idx = i0:i1;

if isempty(bsl_idx)
    warning('Empty baseline index range; no correction applied.');
    EEG.data = X; if ndims(EEG.data)==2, EEG.data = squeeze(EEG.data); end
    return
end

fprintf('Baseline: method="%s", mode="%s", window=[%g %g] ms (idx %d:%d)\n', ...
    method, mode, bsl_ms(1), bsl_ms(2), i0, i1);

% ---------- apply correction ----------
switch method
    case {"median","mean","trimmed mean"}
        % compute baseline statistic per channel × epoch
        switch method
            case "median"
                bsl_vals = median(X(:,bsl_idx,:), 2);   % [C x 1 x N]
            case "mean"
                bsl_vals = mean(X(:,bsl_idx,:), 2);
            case "trimmed mean"
                trim_pct = 20; % %
                trim_amt = round((trim_pct/100) * numel(bsl_idx));
                bsl_vals = nan(C,1,N);
                for ch = 1:C
                    for ep = 1:N
                        vals = squeeze(X(ch,bsl_idx,ep));
                        vals = sort(vals);
                        if numel(vals) > 2*trim_amt
                            vals = vals(trim_amt+1:end-trim_amt);
                        end
                        bsl_vals(ch,1,ep) = mean(vals,'omitnan');
                    end
                end
        end

        switch mode
            case "subtract"
                X = X - repmat(bsl_vals, [1 T 1]);
            case "divide"
                denom = max(eps, abs(bsl_vals));
                X = bsxfun(@rdivide, X, denom);
        end

    case "1/f"
        fs = EEG.srate;
        % frequency vector for n-point FFT on epoch length
        f = (0:T-1)' * (fs/T);
        exclude_band = [4 20]; % Hz to exclude from fit
        fit_mask = f > 0 & (f < exclude_band(1) | f > exclude_band(2)) & f <= fs/2;

        for ch = 1:C
            for ep = 1:N
                x  = squeeze(X(ch,:,ep)).';   % column [T x 1]
                xb = x(bsl_idx);

                Pb = abs(fft(xb, T)).^2;
                if nnz(fit_mask) < 5
                    warning('1/f: too few points for fit (ch=%d, ep=%d). Skipping.', ch, ep);
                    continue;
                end

                % Fit slope/intercept in log-log
                f_log = log10(f(fit_mask));
                P_log = log10(Pb(fit_mask)+eps);
                p_fit = polyfit(f_log, P_log, 1);
                slope = p_fit(1); intercept = p_fit(2);

                % Predict broadband spectrum for full trial
                P1f_pred = 10.^(intercept + slope * log10(f + eps));

                % Subtract broadband magnitude and reconstruct
                X_full = fft(x, T);
                mag = abs(X_full);
                phase = angle(X_full);
                mag_corr = mag - sqrt(P1f_pred(:));
                mag_corr(mag_corr < 0) = 0;

                X_corr = mag_corr .* exp(1i*phase);
                x_corr = real(ifft(X_corr, T));
                X(ch,:,ep) = x_corr;
            end
        end
end

% ---------- write back ----------
EEG.data = X;
end
