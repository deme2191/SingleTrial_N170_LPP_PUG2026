% ============================================================
% find_peak_window.m
% ============================================================
% PURPOSE:
%   Find the peak amplitude and its latency INSIDE a specified
%   time window, using the desired polarity:
%     - 'neg' = most negative point (e.g., N170)
%     - 'pos' = most positive point (e.g., LPP peak)
%
% INPUTS:
%   signal   : [time x 1] vector (a waveform)
%   times_ms : [time x 1] vector of timepoints in milliseconds (EEG.times)
%   win_ms   : [start end] time window in ms, e.g. [140 220]
%   polarity : 'neg' or 'pos'
%
% OUTPUTS:
%   lat_ms   : latency (ms) of the detected peak
%   amp      : amplitude (uV or ICA units) at that peak
% ============================================================
function [lat_ms, amp] = find_peak_window(signal, times_ms, win_ms, polarity)

% Find all sample indices whose time falls inside the requested window
idx = find(times_ms >= win_ms(1) & times_ms <= win_ms(2));

% If the window does not overlap with the epoch time range, return NaNs
if isempty(idx)
    lat_ms = NaN;
    amp    = NaN;
    return;
end

% Extract just the segment of the signal within the time window
seg = signal(idx);

% Depending on polarity, find the minimum (negative peak) or maximum (positive peak)
if strcmpi(polarity,'neg')
    [amp, rel] = min(seg);   % amp = minimum value, rel = index within seg
else
    [amp, rel] = max(seg);   % amp = maximum value, rel = index within seg
end

% Convert the segment-relative index into the original time vector
lat_ms = times_ms(idx(rel));

end
