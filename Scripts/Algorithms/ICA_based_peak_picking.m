% % ICA_based_peak_picking.m
% Authors: Ilayda Cengiz, Zehra Gökce Yilmaz
% Affiliation : Carl von Ossietzky Universität Oldenburg
% Last updated : 2026
%
% PURPOSE:
%   Implements an extended ICA with PCA rank adjustment and single-trial
%   peak-picking to extract N170 and LPP amplitudes and latencies from
%   pre-processed, epoched EEG datasets.
%
%   Expects EEGLAB .set files named according to the fork convention:
%       fork1-2_<ParticipantID>_b<Baseline><Condition>r<Reference>.set
%
%   For each file (participant × condition × baseline × reference):
%     (1) Load the epoched EEGLAB dataset
%     (2) Run a "quantification ICA" (separate from any artifact-removal ICA)
%         with PCA dimensionality reduction matched to the data rank
%     (3) Identify an N170-like IC and an LPP-like IC by correlating
%         component time-courses with sensor-based ROI averages,
%         restricted to the relevant time windows
%     (4) Back-project each selected IC to scalp space and extract:
%             N170 — negative peak amplitude & latency (at P10)
%             LPP  — mean amplitude & positive-peak latency (CP1/CP2/Pz/P3/P4)
%     (5) Save one long-format CSV (one row per trial) per run
%
%   Additional features:
%     - Randomly samples 10% of participants for topography verification
%       and ERP sanity-check figures (saved as .fig and .png)
%     - Labels the N170 and LPP ICs on verification topography plots
%     - Uses ICA convergence tolerance of 1e-3 (optimised for speed)
%     - Each run is saved in a timestamped subfolder to prevent overwrites
%
% DEPENDENCIES:
%   - EEGLAB (tested with eeglab2025.1.0)
%       https://sccn.ucsd.edu/eeglab/
%   - Helper functions (must be on the MATLAB path):
%       find_peak_window.m  — finds a peak (positive or negative) within
%                             a specified time window
%       extract_lpp_metrics.m — extracts LPP mean amplitude and peak latency
%
% USAGE:
%   Edit Section 1 (PATHS) to point to your local EEGLAB installation,
%   helper-function folder, data folder, and output folder, then run.
%
% OUTPUTS:
%   <outPath>/ICA_Results_<timestamp>.csv
%       Long-format table; one row per trial with columns:
%       participant, baseline, reference, condition, data_rank,
%       trial, N170_amp, N170_lat, LPP_amp, LPP_lat
%
%   <outPath>/Verification_Topographies/
%       Topo_<ID>_b<bsl>_<ref>.fig / .png  (sampled participants only)
%
%   <outPath>/Sanity_Checks_ERP/
%       SanityCheck_<ID>_b<bsl>_<ref>.png  (sampled participants only)
%
% REPRODUCIBILITY:
%   A fixed random seed (2026, 'twister') is set at the top of the script
%   so the 10% verification sample is identical across runs on the same
%   participant list.
%
% =============================================================================
clear; clc;
rng('default'); 
rng(2026, 'twister');% Fixed seed for reproducible 10% participant selection

%% 1) PATHS 
%  Edit these four variables to match your local environment.

% Path to your EEGLAB installation
eeglabPath = '/Users/ilaydacengiz/Downloads/eeglab2025.1.0';
 
% Path to the folder containing find_peak_window.m and extract_lpp_metrics.m
helperPath = '/Users/ilaydacengiz/Downloads/Peak-picking_scripts';
 
% Path to the folder containing the pre-processed fork1-2_*.set files
dataPath = '/Users/ilaydacengiz/Downloads/SingleTrial_3participants';
 
% Root output folder; a timestamped subfolder is created automatically
baseOutPath = '/Users/ilaydacengiz/Library/CloudStorage/OneDrive-CarlvonOssietzkyUniversitätOldenburg/Single Trial Project/Peak-picking';
 
% Add dependencies to MATLAB path
addpath(eeglabPath);
addpath(helperPath);
 
% Create a timestamped output directory for this run (avoids overwriting)
runTimestamp = datestr(now, 'yyyymmdd_HHMMSS');
outPath      = fullfile(baseOutPath, ['Run_' runTimestamp]);
topoPath     = fullfile(outPath, 'Verification_Topographies');
sanityPath   = fullfile(outPath, 'Sanity_Checks_ERP');
 
if ~exist(outPath,    'dir'), mkdir(outPath);    end
if ~exist(topoPath,   'dir'), mkdir(topoPath);   end
if ~exist(sanityPath, 'dir'), mkdir(sanityPath); end
 
fprintf('>>> Results for this run will be saved in:\n    %s\n\n', outPath);

%% 2) FIND ALL DATA FILES & SAMPLE PARTICIPANTS 
% Collect all fork1-2_*.set files
files = dir(fullfile(dataPath, 'fork1-2_*.set'));
if isempty(files)
    error('No fork1-2_*.set files found in:\n  %s', dataPath);
end
 
% Parse unique participant IDs from filenames
all_filenames = {files.name};
tokens        = regexp(all_filenames, 'fork1-2_(.*?)_b', 'tokens');
all_subjs     = cellfun(@(x) x{1}{1}, tokens, 'UniformOutput', false);
unique_subjs  = unique(all_subjs);
 
% Randomly sample 10% of participants for topography + sanity-check export.
% ceil() guarantees at least 1 participant is selected even in small datasets.
num_to_sample  = ceil(length(unique_subjs) * 0.10);
sampled_subjs  = unique_subjs(randperm(length(unique_subjs), num_to_sample));
 
fprintf('Total participants : %d\n', length(unique_subjs));
fprintf('Sampled for verification : %d\n', num_to_sample);
fprintf('Sampled IDs : %s\n\n', strjoin(sampled_subjs, ', '));

%% 3)ANALYSIS PARAMETERS 
%  Adjust these to match your experimental design and ERP components of interest.

% N170 time windows (ms) 
% Broad window used to locate the subject-average N170 peak.
N170_search_win = [155 210];
 
% Half-width of the subject-centred trial-level search window.
% Final trial window = [subj_peak - N170_halfwin, subj_peak + N170_halfwin]
N170_halfwin    = 40;           % ms
 
% LPP time window (ms) 
% Single window used for both IC selection and trial-level metric extraction.
LPP_win         = [400 600];
 
% Scalp ROI channel labels
% Used ONLY to identify which IC best represents each component;
% back-projection to individual channels / channel means is done per trial.
N170_chans = {'P10'};                          % right occipito-temporal
LPP_chans  = {'CP1','CP2','Pz','P3','P4'};    % centro-parietal
 
% Representative emotion for verification figures 
% Topography and sanity-check figures are exported only for this condition
% to avoid generating an excessive number of files.
rep_emotion = 'happiness';


%% 4) STORAGE FOR RESULTS 

results = [];  % struct array — one element per trial
row = 1; % running row index

%% 5) MAIN LOOP: ITERATE OVER ALL DATASET FILES
for f = 1:length(files)
    % Load the epoched dataset 
    EEG = pop_loadset('filename',files(f).name,'filepath',dataPath);
    
    % Parse key metadata from the filename 
      %  Expected format:
    %      fork1-2_<ParticipantID>_b<Baseline><Condition>r<Reference>.set
    %
    %  The greedy (.*) match for <Condition> ensures that condition labels
    %  ending in 'r' (e.g. 'anger') are captured fully.
    tok = regexp(files(f).name, 'fork1-2_(.*?)_b(-?\d+)(.*)r(Avg|Mas|CSD)\.set', 'tokens');
    
    subj = tok{1}{1};              % Participant ID (string)          
    bsl  = str2double(tok{1}{2});  % Baseline (ms, numeric)
    cond = tok{1}{3};              % Condition label (string) 
    ref  = tok{1}{4};              % Reference scheme: : Avg | Mas | CSD
    
    fprintf('Processing  %s | cond: %s | baseline: %d ms | ref: %s  [%d / %d]\n', ...
            subj, cond, bsl, ref, f, length(files));

    %  RUN QUANTIFICATION ICA WITH PCA RANK ADJUSTMENT
    %  This ICA is used solely to decompose the EEG into statistically
    %  independent components for single-trial metric extraction.
    %  It is run AFTER any artifact-removal ICA that may have been applied
    %  during preprocessing.
    %
    %  Rank adjustment (via the 'pca' option in runica) prevents runica from
    %  searching for non-existent components. this is essential when data
    %  rank is reduced after interpolation or prior ICA removal.
    %
    %  Convergence tolerance 1e-3 provides a good speed-accuracy trade-off 
    %  for large datasets.
    %  ------------------------------------------

    %  EEG.data is [channels x time x trials]
    %  ICA expects a 2D matrix [channels x (time*trials)]
    %  Reshape 3D data to 2D
    data2d = reshape(EEG.data, size(EEG.data,1), []);

    % Compute effective rank to avoid mathematical singularities
    dataRank = rank(data2d); 
    fprintf('Rank for the current file: %d\n', dataRank);
    
    % Run extended infomax ICA with PCA dimensionality reduction
    [weights, sphere] = runica(data2d, ...
                               'extended', 1,        ... % extended ICA (handles sub-Gaussian sources)
                               'stop',     1e-3,     ... % convergence tolerance
                               'pca',      dataRank);    % reduce to effective rank
 
    % Compute IC activations and reshape back to 3D [ICs × time × trials]
    acts = weights * sphere * data2d;
    acts = reshape(acts, size(acts, 1), EEG.pnts, EEG.trials);

    %  SELECT THE N170-LIKE IC  
    %
    %  Correlate each IC's trial-averaged time-course with the
    %  trial-averaged N170 ROI signal, but only within the N170 search window.
    %  Window-restricted correlation focuses the selection on the component of
    %  interest and reduces spurious matches from other time periods.
    % ------------------------------------

    % Find channel indices for the N170 ROI
    chanIdxN = find(ismember({EEG.chanlocs.labels}, N170_chans));
    
    % Compute trial-averaged ROI signal
    topoN170 = squeeze(mean(mean(EEG.data(chanIdxN,:,:),3),1)); % [time x 1]

    % Define the time index for the N170 window 
    idxN = EEG.times >= N170_search_win(1) & EEG.times <= N170_search_win(2);

   % Correlate each IC's mean timecourse with ROI summary ONLY within the window
    corrN = zeros(size(acts,1),1);
    for ic = 1:size(acts,1)
        ic_tc = squeeze(mean(acts(ic,:,:),3)); % mean across trials
        % Only correlate the samples inside the window
        corrN(ic) = corr(ic_tc(idxN)', topoN170(idxN)', 'rows','complete');
    end

    % Choose the IC with maximum correlation
    [~, N170_ic] = max(corrN);

    % SELECT THE LPP-LIKE IC 
    % Same strategy, but using an LPP ROI and the LPP window
    chanIdxL = find(ismember({EEG.chanlocs.labels}, LPP_chans));
    topoLPP = squeeze(mean(mean(EEG.data(chanIdxL,:,:),3),1)); % [time x 1]

    % Define the time index for the LPP window 
    idxL = EEG.times >= LPP_win(1) & EEG.times <= LPP_win(2);

    corrL = zeros(size(acts,1),1);
    for ic = 1:size(acts,1)
        ic_tc = squeeze(mean(acts(ic,:,:),3)); % mean across trials
        % Only correlate the samples inside the LPP window
        corrL(ic) = corr(ic_tc(idxL)', topoLPP(idxL)', 'rows','complete');
    end
    [~, LPP_ic] = max(corrL);
    % AUTOMATED VERIFICATION EXPORT (.fig & .png) -----
    % Logic: Generate topographies only for sampled participants and the representative condition
    if ismember(subj, sampled_subjs) && strcmpi(cond, rep_emotion)
        
        % Construct filenames for both .fig and .png formats
        figName = sprintf('Topo_%s_b%d_%s', subj, bsl, ref);
        fullFigPath = fullfile(topoPath, [figName '.fig']);
        fullPngPath = fullfile(topoPath, [figName '.png']); 
        
        % Proceed only if the figure hasn't been generated yet (to save time)
        if ~exist(fullFigPath, 'file')
            
            % Calculate the mixing matrix (Topographies) from ICA weights
            topo_map = pinv(weights * sphere);
            
            % Create figure (Visible off keeps the UI clean during the long run)
            h = figure('Name', figName, 'Visible', 'off', 'Color', 'w'); 
            
            % --- GLOBAL FIGURE TITLE ---
            % Includes Participant ID, Baseline, and Reference
            % 'Interpreter', 'none' prevents underscores from being rendered as italics
            sgtitle(sprintf('Participant: %s | Baseline: %d ms | Reference: %s', ...
                    subj, bsl, ref), 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
            
            % Plot the top 20 Independent Components
            for j = 1:20
                subplot(5, 4, j); 
                topoplot(topo_map(:,j), EEG.chanlocs);
                
                % IC Labeling: Identify which IC was chosen for N170 and LPP
                title_str = ['IC ' num2str(j)];
                if j == N170_ic, title_str = [title_str ' [N170]']; end
                if j == LPP_ic,  title_str = [title_str ' [LPP]']; end
                
                title(title_str, 'FontSize', 9, 'FontWeight', 'normal');
            end
            
            % Supplemental: Plot N170 if it is not in the first 20
            if N170_ic > 20
                subplot(6, 4, 21);
                topoplot(topo_map(:, N170_ic), EEG.chanlocs);
                title(['IC ' num2str(N170_ic) ' [N170]'], 'FontSize', 9, 'Color', 'r');
            end
            
            % Supplemental: Plot LPP if it is not in the first 20
            if LPP_ic > 20
                subplot(6, 4, 22);
                topoplot(topo_map(:, LPP_ic), EEG.chanlocs);
                title(['IC ' num2str(LPP_ic) ' [LPP]'], 'FontSize', 9, 'Color', 'r');
            end
            
            set(h, 'Visible', 'on');

            % --- SAVE SECTION ---
            % Save as MATLAB .fig for future editing
            saveas(h, fullFigPath); 
            % Save as .png for quick review, reports, and presentations
            saveas(h, fullPngPath); 
            
            % Close figure to free up system memory
            close(h); 
            fprintf('>>> Verification Exported (FIG & PNG): %s\n', figName);
        end
    end
    % DEFINE SUBJECT-SPECIFIC N170 WINDOW 
    %
    % For N170, peak-picking is more reliable if we center the window
    % around each subject's average N170 latency.
    subj_avg_N = squeeze(mean(acts(N170_ic,:,:),3));  % average across trials => [time x 1]

    % Find the subject's N170 peak latency in a broad window
    [lat0, ~] = find_peak_window(subj_avg_N, EEG.times(:), N170_search_win, 'neg');

    % Build a narrower trial-level window centered on subject latency
    N170_win = [lat0 - N170_halfwin, lat0 + N170_halfwin];

    % SINGLE TRIAL METRIC EXTRACTION
    %
    %  Both components are back-projected from IC space to scalp space using
    %  the mixing matrix before peak-picking. This preserves the original
    %  amplitude scale (µV or µV/m²) at the scalp sensors.
    %
    %  N170: negative peak amplitude and latency at P10
    %  LPP : mean amplitude and positive-peak latency over the centro-parietal ROI
    % -------------------------------------------------------------------------
    
    % Calculate mixing matrix (pseudo-inverse of weights*sphere)
    mixingMatrix = pinv(weights * sphere); 

    % Channel indices for back-projected scalp signals
    idx_P10 = find(strcmpi('P10', {EEG.chanlocs.labels})); 
    idx_LPP_roi = find(ismember({EEG.chanlocs.labels}, LPP_chans));

    for tr = 1:EEG.trials
        % N170: Full Back-projection to Scalp, then pick P10 
        sigN_ic = squeeze(acts(N170_ic,:,tr));  % [1 × time]
        scalp_recon_N = mixingMatrix(:,N170_ic) * sigN_ic; % [channels × time]
        sigN_uv = scalp_recon_N(idx_P10,:)'; % [time × 1], P10 only% [time × 1], P10 only
        
        [n170_lat, n170_amp] = find_peak_window(sigN_uv, EEG.times(:), N170_win, 'neg');
        
        % LPP: Full Back-projection to Scalp, then average ROI 
        sigL_ic = squeeze(acts(LPP_ic,:,tr)); % [1 × time]
        scalp_recon_L = mixingMatrix(:,LPP_ic) * sigL_ic; % [channels × time]
        scalp_roi = scalp_recon_L(idx_LPP_roi,:); % [ROI chans × time]
        sigL_uv = mean(scalp_roi,1)';    % [time × 1], ROI mean

        [lpp_amp, lpp_lat] = extract_lpp_metrics(sigL_uv, EEG.times(:), LPP_win);

        % Store results in struct 
        results(row).participant = subj;
        results(row).baseline    = bsl;
        results(row).reference   = ref;
        results(row).condition   = cond;
        results(row).data_rank   = dataRank;
        results(row).trial       = tr;

        results(row).N170_amp = n170_amp;
        results(row).N170_lat = n170_lat;

        results(row).LPP_amp  = lpp_amp;
        results(row).LPP_lat  = lpp_lat;

        row = row + 1;
    end

    %  ERP SANITY CHECKS (sampled participants × representative emotion)
    %
    %  Validates that the full reconstruction from all ICs (A × activations)
    %  closely matches the original ERP at key channels / ROIs.
    %  A large discrepancy would indicate a problem with the ICA solution.
    %
    %  Subplot 1  N170: original vs. reconstructed ERP at P10
    %  Subplot 2  LPP : original vs. reconstructed ERP (ROI mean)
    % -------------------------------------------------------------------------
   
    if ismember(subj, sampled_subjs) && strcmpi(cond, rep_emotion)
        
        sanityFigName = sprintf('SanityCheck_%s_b%d_%s', subj, bsl, ref);
        fullSanityPath = fullfile(sanityPath, [sanityFigName '.png']);
        
        if ~exist(fullSanityPath, 'file')
            % Full reconstruction: apply mixing matrix to all IC activations
            Xrec_2d = mixingMatrix * reshape(acts, size(acts,1), []);
            Xrec_3d = reshape(Xrec_2d, size(EEG.data));
            
            % Trial-averaged ERPs
            orig_erp  = mean(EEG.data, 3);
            recon_erp = mean(Xrec_3d, 3);
            
            % Y-axis label reflects the reference-dependent amplitude unit
            if strcmpi(ref, 'CSD')
                y_label_str = '\muV / normalised surface area^2';
            else
                y_label_str = '\muV';
            end
           
            h_san = figure('Name', sanityFigName, 'Visible', 'off', 'Color', 'w');
            
            % Subplot 1: N170 (P10 Single Channel) 
            subplot(2, 1, 1);
            idx_p10 = find(strcmpi('P10', {EEG.chanlocs.labels}));
            plot(EEG.times, orig_erp(idx_p10, :), 'k'); hold on;
            plot(EEG.times, recon_erp(idx_p10, :), 'r--');
            title(['N170 Check: P10 Channel']);
            ylabel('\muV'); legend('Original', 'Reconstructed'); grid on;

            % Subplot 2: LPP (ROI AVERAGE) 
            subplot(2, 1, 2);
            % Using the ROI channels we specified (CP1, CP2, Pz, P3, P4)
            idx_lpp_roi = find(ismember({EEG.chanlocs.labels}, LPP_chans));
            
            % Calculate Mean of the ROI for both original and reconstructed ERP
            orig_lpp_roi_mean = mean(orig_erp(idx_lpp_roi, :), 1);
            recon_lpp_roi_mean = mean(recon_erp(idx_lpp_roi, :), 1);
            
            plot(EEG.times, orig_lpp_roi_mean, 'k'); hold on;
            plot(EEG.times, recon_lpp_roi_mean, 'r--');
            title(['LPP Check: ROI Mean (CP1, CP2, Pz, P3, P4)']);
            xlabel('Time (ms)'); ylabel('\muV'); grid off;
            
            sgtitle(['Sanity Check: ' subj ' (' ref ')'], 'Interpreter', 'none');
            
            saveas(h_san, fullSanityPath);
            close(h_san);
            fprintf('>>> Sanity Check (ROI-based) Saved: %s\n', sanityFigName);
        end
    end
end

%% 6) SAVE RESULTS 
% Convert struct array to table and write CSV
T = struct2table(results);
outFile = fullfile(outPath, sprintf('ICA_Results_%s.csv', runTimestamp));
writetable(T, outFile, 'Delimiter', ';');

fprintf('ICA processing complete.\n');
fprintf('CSV Results: %s\n', outFile);
fprintf('Verification Figures: %s\n', topoPath);
fprintf('  Topographies       : %s\n', topoPath);
fprintf('  ERP sanity checks  : %s\n', sanityPath);
