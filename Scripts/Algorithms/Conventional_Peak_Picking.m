% ============================================================
% Conventional_Peak_Picking.m
%
% Authors : Ilayda Cengiz, Zehra Gökce Yilmaz
% Affiliation : Carl von Ossietzky Universität Oldenburg
% Last updated : 2026
%
% PURPOSE:
%   Implements conventional (ICA-free) single-trial peak-picking directly on
%   scalp electrode data. 
%
%   Expects EEGLAB .set files named according to the fork convention:
%       fork1-2_<ParticipantID>_b<Baseline><Condition>r<Reference>.set
%
%   For each file (participant × condition × baseline × reference):
%     (1) Load the epoched EEGLAB dataset
%     (2) Average data across predefined ROI channels to obtain one
%         waveform per trial (no source separation performed)
%     (3) Locate each participant's average N170 peak to define a
%         subject-centred trial-level search window 
%     (4) Extract single-trial metrics:
%             N170 — negative peak amplitude & latency (P10)
%             LPP  — mean amplitude & positive-peak latency (CP1/CP2/Pz/P3/P4)
%     (5) Save one long-format CSV (one row per trial) per run
%
% DEPENDENCIES:
%   - Helper functions (must be on the MATLAB path):
%       find_peak_window.m    — finds a peak (positive or negative) within
%                               a specified time window
%       extract_lpp_metrics.m — extracts LPP mean amplitude and peak latency
%
% OUTPUTS:
%   <outPath>/Simple_Scalp_Peak_Picking_<timestamp>.csv
%       Long-format table; one row per trial with columns:
%       participant, baseline, reference, condition, algorithm,
%       trial, N170_amp, N170_lat, LPP_amp, LPP_lat
% ====================================================

clear; clc;

%% 1) PATHS 
%  Edit these four variables to match your local environment.

% Path to your EEGLAB 
eeglabPath = '/Users/ilaydacengiz/Downloads/eeglab2025.1.0';
 
% Path to the folder containing find_peak_window.m and extract_lpp_metrics.m
helperPath = '/Users/ilaydacengiz/Downloads/Peak-picking_scripts';
 
% Path to the folder containing the preprocessed fork1-2_*.set files
dataPath = '/Users/ilaydacengiz/Downloads/erp';

% Add dependencies to MATLAB path
addpath(eeglabPath);
addpath(helperPath);
 
% Output folder for the results CSV
outPath = '/Users/ilaydacengiz/Library/CloudStorage/OneDrive-CarlvonOssietzkyUniversitätOldenburg/Single Trial Project/Peak-picking';
if ~exist(outPath,'dir'), mkdir(outPath); end


%% 2) FIND FILES 
files = dir(fullfile(dataPath, 'fork1-2_*.set'));
 
if isempty(files)
    error('No fork1-2_*.set files found in:\n  %s', dataPath);
end
 
fprintf('Total files found: %d. Starting analysis...\n\n', length(files));

%% 3) PARAMETERS 
% N170: occipito-temporal
N170_search_win = [155 210];   

% Half-width of the subject-centred trial-level search window.
% Final trial window = [subj_peak - N170_halfwin, subj_peak + N170_halfwin]
N170_halfwin    = 40;          
N170_chans      = {'P10'};      

% LPP: centro-parietal 
LPP_win         = [400 600];
LPP_chans       = {'CP1','CP2','Pz','P3','P4'};       

%% 4) STORAGE
results = [];
row = 1;

%% 5) MAIN LOOP 
for f = 1:length(files)
    EEG = pop_loadset('filename', files(f).name, 'filepath', dataPath);
    
    % Parse metadata
    tok = regexp(files(f).name, 'fork1-2_(.*?)_b(-?\d+)(.*)r(Avg|Mas|CSD)\.set', 'tokens');
    subj = tok{1}{1}; bsl = str2double(tok{1}{2}); cond = tok{1}{3}; ref = tok{1}{4};
    
    fprintf('Processing %s | %s | b%d | %s (File %d/%d)\n', subj, cond, bsl, ref, f, length(files));

    %% PREPARE ROI DATA
    %
    %  Unlike the ICA-based script, no source separation is applied.
    %  Instead, signals are averaged across the predefined ROI electrodes to
    %  yield a single representative waveform per trial per component.
    %
    %  EEG.data is [channels × time × trials].
    %  After averaging across channels and squeezing, roi_data* is [time × trials].
    % --------------------------------------
 
    chanIdxN = find(ismember({EEG.chanlocs.labels}, N170_chans));
    chanIdxL = find(ismember({EEG.chanlocs.labels}, LPP_chans));
 
    % Average over ROI channels → [time × trials]
    % squeeze() removes the singleton channel dimension after mean()
    roi_dataN = squeeze(mean(EEG.data(chanIdxN, :, :), 1));
    roi_dataL = squeeze(mean(EEG.data(chanIdxL, :, :), 1));

    %% SUBJECT-SPECIFIC N170 WINDOW 
     % Trial-average of the N170 ROI waveform → [time × 1]
    subj_avg_N = mean(roi_dataN, 2);
 
    % Locate the subject's N170 peak within the broad search window
    [lat0, ~] = find_peak_window(subj_avg_N, EEG.times(:), N170_search_win, 'neg');
 
    % Subject-centred trial-level search window
    N170_win = [lat0 - N170_halfwin, lat0 + N170_halfwin];
    %% SINGLE-TRIAL METRIC EXTRACTION
       for tr = 1:EEG.trials
 
        % N170: negative peak amplitude and latency at P10 (or P10 ROI mean)
        sigN = roi_dataN(:, tr);    % [time × 1]
        [n170_lat, n170_amp] = find_peak_window(sigN, EEG.times(:), N170_win, 'neg');
 
        % LPP: mean amplitude and positive-peak latency over centro-parietal ROI
        sigL = roi_dataL(:, tr);    % [time × 1]
        [lpp_amp, lpp_lat] = extract_lpp_metrics(sigL, EEG.times(:), LPP_win);
 
        % Store result (one row per trial)
        results(row).participant = subj;
        results(row).baseline    = bsl;
        results(row).reference   = ref;
        results(row).condition   = cond;
        results(row).algorithm   = 'Simple_Scalp_Peak';  % identifier for cross-method comparison
        results(row).trial       = tr;
        results(row).N170_amp    = n170_amp;
        results(row).N170_lat    = n170_lat;
        results(row).LPP_amp     = lpp_amp;
        results(row).LPP_lat     = lpp_lat;
 
        row = row + 1;
    end
 
end % end main file loop

%% SAVE RESULTS
% Convert struct array to table
T = struct2table(results);
 
outFile = fullfile(outPath, ...
          sprintf('Simple_Scalp_Peak_Picking_%s.csv', datestr(now, 'yyyymmdd_HHMMSS')));
writetable(T, outFile, 'Delimiter', ';');
 
fprintf('\n============================================================\n');
fprintf('Processing complete.\n');
fprintf('  CSV results : %s\n', outFile);
