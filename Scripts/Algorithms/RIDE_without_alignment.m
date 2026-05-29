%% RIDE_without_alignment.m
%
% Authors : Ilayda Cengiz, Zehra Gökce Yilmaz
% Affiliation : Carl von Ossietzky Universität Oldenburg
% Last updated : 2026
%
% PURPOSE:
%   Implements RIDE (Residue Iteration Decomposition)-based single-trial
%   peak-picking for N170 and LPP.nRIDE decomposes the EEG into latency-variable 
%   components without assuming fixed temporal structure, making it 
%   well-suited for single-trial analysis of components with trial-to-trial 
%   latency jitter.
%
%     - N170 is extracted from the S-component (stimulus-locked)
%     - LPP is extracted from the C-component (response-locked / variable)
%       WITHOUT synchronisation, so natural trial-to-trial latency variance
%       is preserved.
%   Expects EEGLAB .set files named according to the fork convention:
%       fork1-2_<ParticipantID>_b<Baseline><Condition>r<Reference>.set
%
%   For each file (participant × condition × baseline × reference):
%     (1) Load the epoched EEGLAB dataset
%     (2) Permute data to RIDE's expected format [time × channels × trials]
%     (3) Run RIDE decomposition (S and C components)
%     (4) Isolate each component at the relevant ROI channels using
%         single_trial_RIDE.m and average across ROI channels
%     (5) Extract single-trial metrics:
%             N170 — negative peak amplitude & latency from S-component (P10)
%             LPP  — mean amplitude & positive-peak latency from C-component
%                    (CP1/CP2/Pz/P3/P4), unsynchronised
%     (6) Save one long-format CSV (one row per trial) per run
% DEPENDENCIES:
%   - EEGLAB (tested with eeglab2025.1.0)
%       https://sccn.ucsd.edu/eeglab/
%   - RIDE toolbox (RIDE-main)
%       https://github.com/guangouyang/RIDE  (or your local copy)
%       Required functions: RIDE_cfg.m, RIDE_call.m, single_trial_RIDE.m
%   - Helper functions (must be on the MATLAB path):
%       find_peak_window.m    — finds a peak (positive or negative) within
%                               a specified time window
%       extract_lpp_metrics.m — extracts LPP mean amplitude and peak latency
%
% USAGE:
%   Edit Section 1 (PATHS) to point to your local toolbox installations,
%   data folder, and output file path, then run.
%
% OUTPUTS:
%   RIDE_Results_<timestamp>.csv  (written to the working directory or outPath)
%       Long-format table; one row per trial with columns:
%       participant, baseline, reference, condition, trial,
%       N170_amp, N170_lat, LPP_amp, LPP_lat
% ================================================


clear; clc;

%% PATHS AND CONFIGURATION
%  Edit these variables to match your local environment.

% Path to the folder containing the pre-processed fork1-2_*.set files
input_dir = 'C:\Users\ilaydaandzehra\Downloads\Ilayda-and-Zehra';
 
% Path to your RIDE toolbox (genpath adds all subdirectories automatically)
ridePath = 'C:\Users\ilaydaandzehra\Downloads\RIDE-main';
 
% Path to your EEGLAB installation
eeglabPath = 'C:\Users\ilaydaandzehra\Downloads\eeglab2025.1.0';
 
% Path to the folder containing find_peak_window.m and extract_lpp_metrics.m
helperPath = 'C:\Users\ilaydaandzehra\Downloads\Peak-picking_scripts';
 
% Output CSV path (written to the current working directory)
outFile = ['RIDE_Results_', datestr(now, 'yyyymmdd_HHMMSS'), '.csv'];
 
% Add all dependencies to the MATLAB path
addpath(genpath(ridePath));
addpath(eeglabPath);
addpath(helperPath);
 
eeglab; close all;

%% Analysis Parameters
% N170 time window (ms) 
N170_win   = [155 210]; 

% LPP time window (ms)
LPP_win    = [400 600];

% Scalp ROI channel labels
N170_chans = {'P10'};   
LPP_chans  = {'CP1','CP2','Pz','P3','P4'}; 

% Condition filter
% This script currently processes only the happiness condition.
% Change or remove the filter in the next section to process all conditions.
target_condition = 'happiness';

% Find all data files
files = dir(fullfile(input_dir, 'fork1-2_*.set'));

if isempty(files)
    error('No fork1-2_*.set files found in:\n  %s', input_dir);
end

% Inıtıalise results storage
results_all = [];  % struct array — one element per trial
row_idx = 1; % running row index

%% MAIN LOOP
for f = 1:length(files)
    % Condition filter: Skip files that do not match the target condition. 
    if ~contains(lower(files(f).name), lower(target_condition))
        continue;
    end
    
    fprintf('Processing: %s\n', files(f).name);

    % Load dataset
    EEG = pop_loadset('filename', files(f).name, 'filepath', char(input_dir));
    
    % Metadata Extraction
    %  Expected format:
    %      fork1-2_<ParticipantID>_b<Baseline><Condition>r<Reference>.set
    tok = regexp(files(f).name, 'fork1-2_(.*?)_b(-?\d+)(.*?)r(Avg|Mas|CSD)\.set', 'tokens');
    if isempty(tok)
        warning('Could not parse filename: %s — skipping.', files(f).name);
        continue;
    end

    subj = tok{1}{1};   % Participant ID  (string)
    bsl  = tok{1}{2};   % Baseline onset  (string; convert if numeric needed)
    cond = tok{1}{3};   % Condition label (string)
    ref  = tok{1}{4};   % Reference scheme: Avg | Mas | CSD
 
    % Permute data for RIDE compatibility 
    %  EEGLAB stores data as [channels × time × trials].
    %  RIDE expects [time × channels × trials].
    %  permute() reorders the dimensions without copying data values.
    data_ride = permute(EEG.data, [2 1 3]); % [time × channels × trials]
    t_axis = EEG.times;
    Fs = EEG.srate;
    
    %% RIDE DECOMPOSITION
    %  Two components are requested:
    %    's' (stimulus-locked): used to isolate the N170
    %    'c' (variable/response-locked): used to isolate the LPP
    %
    %  Component time windows (cfg.comp.twd) define the plausible latency
    %  range for each component; latency 0 pins 's' to stimulus onset.
    %
    %  RIDE_cfg() fills in all remaining defaults; RIDE_call() runs the
    %  iterative decomposition.
    % ------------------------------------
    cfg = [];
    cfg.samp_interval = 1000/Fs; % sampling interval (ms)
    cfg.epoch_twd     = [EEG.xmin*1000 EEG.xmax*1000]; % epoch window (ms)
    cfg.comp.name     = {'s', 'c'}; 
    cfg.comp.twd      = {[0 400], [200 800]}; % plausible latency ranges (ms)
    cfg.comp.latency  = {0, 'unknown'}; % 's' is stimulus-locked; 'c' is estimated

    cfg = RIDE_cfg(cfg); % fill in RIDE defaults
    res = RIDE_call(data_ride, cfg);  % run decomposition

    % ROI channel indices 
    idxN = find(ismember({EEG.chanlocs.labels}, N170_chans));
    idxL = find(ismember({EEG.chanlocs.labels}, LPP_chans));

    %% COMPONENT ISOLATION 
    %  single_trial_RIDE() returns the isolated component waveform for one
    %  channel at a time [time × trials]. We loop over ROI channels and
    %  average, consistent with how the other two scripts handle multi-channel
    %  ROIs. For N170_chans = {'P10'} this loop runs once.
    % --------------------------------------
   
    % N170: Isolate 's' (Stimulus-locked)
    st_s_signals = zeros(size(data_ride,1), length(idxN), EEG.trials);
   
    for ch = 1:length(idxN)
        % single_trial_RIDE extracts isolated component based on channel index
        st_s_signals(:,ch,:) = single_trial_RIDE(data_ride, res, 's', idxN(ch));
    end

    % Average over ROI channels :[time × trials]
    sig_N_clean = squeeze(mean(st_s_signals, 2)); % Average if multiple ROI chans

    % LPP Extraction: Isolate 'c' WITHOUT Synchronization 
    %  The C-component is extracted WITHOUT synchronisation so that
    %  natural trial-to-trial latency variance is preserved. Synchronising
    %  would align all trials to the component's estimated latency, which
    %  would artificially inflate peak amplitudes and remove the latency
    %  variability we want to capture.

    st_c_signals = zeros(size(data_ride,1), length(idxL), EEG.trials);
    
    for ch = 1:length(idxL)
        % We remove 'synced' to keep the original trial-by-trial timing
        st_c_signals(:,ch,:) = single_trial_RIDE(data_ride, res, 'c', idxL(ch)); 
    end

     % Average over ROI channels
    sig_L_clean = squeeze(mean(st_c_signals, 2)); % [Time x Trial] - Original timing, but S-free

    %% SINGLE-TRIAL METRIC EXTRACTION
    for tr = 1:EEG.trials
        % Measure N170 from isolated S-component
        [n170_lat, n170_amp] = find_peak_window(sig_N_clean(:, tr), t_axis(:), N170_win, 'neg');

        % Measure LPP from isolated & aligned C-component
        [lpp_amp, lpp_lat] = extract_lpp_metrics(sig_L_clean(:, tr), t_axis(:), LPP_win);

        % Store results
        results_all(row_idx).participant = subj;
        results_all(row_idx).baseline    = bsl;
        results_all(row_idx).reference   = ref;
        results_all(row_idx).condition   = cond;
        results_all(row_idx).trial       = tr;
        results_all(row_idx).N170_amp    = n170_amp;
        results_all(row_idx).N170_lat    = n170_lat;
        results_all(row_idx).LPP_amp     = lpp_amp;
        results_all(row_idx).LPP_lat     = lpp_lat;
        row_idx = row_idx + 1;
    end
end

%% SAVE

if isempty(results_all)
    warning('No results were collected. Check that files matching "%s" exist.', ...
            target_condition);
else
    % Convert struct array to table 
    writetable(struct2table(results_all), outFile, 'Delimiter', ';');
 
    fprintf('Extraction complete.\n');
    fprintf('  CSV results : %s\n', outFile);
    fprintf('============================================================\n');
end