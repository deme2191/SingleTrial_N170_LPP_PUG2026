% ============================================================
% extract_lpp_metrics.m
% ============================================================
% PURPOSE:
%   Extract LPP metrics from a waveform in a specified time window:
%     1) Mean amplitude in the LPP window (recommended for broad components)
%     2) Peak latency (time of maximum positivity) in the same window (optional)
%
% INPUTS:
%   signal   : [time x 1] vector (a waveform)
%   times_ms : [time x 1] timepoints in ms (EEG.times)
%   win_ms   : [start end] time window in ms (e.g., [300 800])
%
% OUTPUTS:
%   mean_amp : mean amplitude over the window
%   peak_lat : latency (ms) of the maximum value within the window
% ============================================================
function [mean_amp, peak_lat] = extract_lpp_metrics(signal, times_ms, win_ms)

% Find indices corresponding to the requested time window
idx = find(times_ms >= win_ms(1) & times_ms <= win_ms(2));

% If no indices fall in the window, return NaNs
if isempty(idx)
    mean_amp = NaN;
    peak_lat = NaN;
    return;
end

% Extract windowed segment
seg = signal(idx);

% Mean amplitude in the window
mean_amp = mean(seg);

% Peak latency = timepoint of the maximum value in the window
[~, rel] = max(seg);
peak_lat = times_ms(idx(rel));

end
