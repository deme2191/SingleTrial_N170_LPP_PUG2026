% =============================================================================
% Variance_Decomposition_Happiness.m
%
% Authors : Ilayda Cengiz, Zehra Gökce Yilmaz
% Affiliation : Carl von Ossietzky Universität Oldenburg
% Last updated : 2026
%
% PURPOSE:
%   Quantifies single-trial reliability of N170 and LPP metrics produced by
%   three peak-picking algorithms (ICA-based, RIDE, Conventional) using a
%   Generalizability Theory (G-theory) framework implemented via Linear Mixed
%   Effects (LME) models.
%
%   For each combination of algorithm × metric × baseline × reference:
%     (1) Load the algorithm's CSV output and filter for the target condition
%     (2) Fit a one-way random-effects LME:
%             metric ~ 1 + (1 | participant)
%         to partition variance into between-person (signal) and
%         within-person (noise / trial-to-trial variability) components
%     (3) Compute:
%             ICC          — single-trial intraclass correlation coefficient
%                            (signal variance / total variance)
%             G-coefficient — person-level dependability accounting for
%                            each participant's actual trial count
%             95% bootstrap CI for ICC (participant-level resampling)
%     (4) Append a summary row to the results table
%   Results are saved to a single semicolon-delimited CSV.
%
% STATISTICAL NOTES:
%   - Only one condition (happiness) is analysed here, so 'condition' is
%     excluded from the LME formula. For multi-condition designs, add
%     '+ condition' as a fixed effect.
%   - Bootstrap resampling is at the participant level (cluster bootstrap)
%     to respect the nested structure of trials within participants.
%     Each bootstrap replicate assigns a unique Pseudo-ID so the LME treats
%     resampled duplicates as independent participants.
%   - If an LME bootstrap replicate fails to converge, that replicate is
%     set to NaN and excluded from the CI percentiles.
%
%
% USAGE:
%   Edit Section 1 (CONFIGURATION) to update file paths and parameters,
%   then run. No other changes should be needed for a standard analysis.
%
% OUTPUTS:
%   <outFile>  (semicolon-delimited CSV)
%       One row per algorithm × metric × baseline × reference combination.
%       Columns: Baseline, Reference, Algorithm, Measurement,
%                Var_Between, Var_Within,
%                Var_Between_Pct, Var_Within_Pct,
%                ICC_SingleTrial, CI_95_Low, CI_95_High,
%                Dependability_G, Avg_Trials
%
% ==============================================

close all; clear all; clc;

%  SECTION 1 — CONFIGURATION
%  Edit file paths and analysis settings here.
% ============================================
 
% Input CSV files (one per algorithm, in the same order as 'algorithms')
% Delimiter is detected automatically per file (see Section 3).
filePaths = { ...
    '/Users/ilaydacengiz/Nextcloud/SingleTrial_PracticalProject/LME Reliability/ICA_Results_20260227_210123.csv', ...
    '/Users/ilaydacengiz/Nextcloud/SingleTrial_PracticalProject/LME Reliability/RIDE_Happiness_Results_with_alignment2026-03-09_1926.csv', ...
    '/Users/ilaydacengiz/Nextcloud/SingleTrial_PracticalProject/LME Reliability/Scalp_Results_Final_AllData_20260212_1458.csv' ...
};
 
% Algorithm label for each file (must match filePaths order)
algorithms = {'ICA', 'RIDE', 'Conventional'};
 
% ERP metrics to analyse
metrics = {'N170_amp', 'N170_lat', 'LPP_amp', 'LPP_lat'};
 
% Baseline conditions and reference schemes to loop over
baselines  = [-100, -200];
references = {'Avg', 'Mas', 'CSD'};
 
% Condition to analyse (lowercase; comparison is case-insensitive)
target_condition = 'happiness';
 
% Number of bootstrap replicates for ICC confidence intervals
nBoot = 1000;
 
% Output CSV file path
outFile = ['Reliability_Report_', datestr(now, 'yyyymmdd_HHMMSS'), '.csv'];
 
%  SECTION 2 — INITIALISE RESULTS TABLE
% ==========================================
 
allResults = table();
fprintf('Initiating LME-based Reliability Analysis...\n');
 
%  SECTION 3 — MAIN LOOP: ALGORITHM × BASELINE × REFERENCE × METRIC
% ==========================================
 
for f = 1:length(filePaths)
 
    currentAlgo = algorithms{f};
    fprintf('\n--- Algorithm: %s ---\n', currentAlgo);
 
    %  3a) LOAD CSV
    % All metric columns are configuredto fill missing values with NaN rather than erroring.

    opts           = detectImportOptions(filePaths{f});
    opts.Delimiter = ';';
    opts           = setvaropts(opts, metrics, 'FillValue', NaN);
    data           = readtable(filePaths{f}, opts);
 
    % Standardise variable types used in filtering
    data.participant = categorical(data.participant);
    data.condition   = lower(strtrim(string(data.condition))); % normalise case/whitespace
 
    %  3b) LOOP OVER BASELINE × REFERENCE COMBINATIONS
    for bl = baselines
        for ref_idx = 1:length(references)
 
            currentRef = references{ref_idx};
 
            % Subset to current baseline and reference
            subData = data(data.baseline == bl & strcmpi(data.reference, currentRef), :);
            if isempty(subData), continue; end
 
            %  3c) LOOP OVER METRICS
            for m = 1:length(metrics)
 
                targetMetric = metrics{m};
                fprintf('  Analyzing  %-12s | BL: %4d ms | Ref: %s\n', ...
                        targetMetric, bl, currentRef);
 
                try
 
                    %  STEP 1 — FILTER: target condition + non-missing metric
                    metricData = subData( ...
                        strcmpi(subData.condition, target_condition) & ...
                        ~isnan(subData.(targetMetric)), :);
 
                    if isempty(metricData)
                        fprintf('    !! Skipping: no %s data for condition "%s"\n', ...
                                targetMetric, target_condition);
                        continue;
                    end
 
                    % Per-participant trial counts (used for G-coefficient)
                    [uniqueParts, ~, idx] = unique(metricData.participant);
                    countsPerPerson       = accumarray(idx, 1);
                    nParts                = numel(uniqueParts);
 
                    %  STEP 2 — FIT LME AND PARTITION VARIANCE
                    %
                    %  Formula: metric ~ 1 + (1 | participant)
                    %    Fixed effect  : grand mean (intercept only)
                    %    Random effect : participant-level deviation
                    %
                    %  covarianceParameters returns:
                    %    psi{1}  — between-person variance (tau²)
                    %    sigma   — residual / within-person variance (sigma²)

                    formula = [targetMetric ' ~ 1 + (1 | participant)'];
                    lme     = fitlme(metricData, formula);
 
                    [psi, sigma] = covarianceParameters(lme);
                    var_between  = double(psi{1});
                    var_within   = double(sigma);
                    total_var    = var_between + var_within;
 
                    %  STEP 3 — COMPUTE RELIABILITY INDICES
                    %
                    %  ICC (single-trial):
                    %    Proportion of total variance attributable to stable
                    %    between-person differences. Values near 1 indicate
                    %    that single trials reliably rank participants.
                    %
                    %  G-coefficient (dependability):
                    %    Person-level reliability accounting for each
                    %    participant's actual trial count. Equivalent to the
                    %    Spearman-Brown corrected reliability for unequal n.
                   
                    icc_single      = var_between / total_var;
                    dep_per_person  = var_between ./ ...
                                      (var_between + (var_within ./ countsPerPerson));
                    dependability_g = mean(dep_per_person);
 
                    %  STEP 4 — BOOTSTRAP 95% CI FOR ICC
                    %
                    %  Participant-level (cluster) bootstrap:
                    %    Resample whole participants with replacement so that
                    %    the nested trial structure is preserved. Each replicate
                    %    assigns unique Pseudo-IDs to avoid duplicate-label
                    %    issues in fitlme.
 
                    % Pre-split data into per-participant cells for speed
                    partCells = cell(nParts, 1);
                    for p = 1:nParts
                        partCells{p} = metricData(metricData.participant == uniqueParts(p), :);
                    end
 
                    bootICCs = NaN(nBoot, 1);
 
                    for b = 1:nBoot
 
                        % Resample participants with replacement
                        resampIdx = randi(nParts, nParts, 1);
                        bootTmp   = cell(nParts, 1);
 
                        for ri = 1:nParts
                            pData = partCells{resampIdx(ri)};
                            % Assign unique Pseudo-ID so LME treats each
                            % resampled copy as a distinct participant
                            pData.participant(:) = categorical({['Pseudo' num2str(ri)]});
                            bootTmp{ri} = pData;
                        end
 
                        bootData = vertcat(bootTmp{:});
 
                        try
                            lme_b            = fitlme(bootData, formula);
                            [psi_b, sigma_b] = covarianceParameters(lme_b);
                            v_b              = max(0, double(psi_b{1})); % clamp: variance >= 0
                            v_w              = double(sigma_b);
                            bootICCs(b)      = v_b / (v_b + v_w);
                        catch
                            % Convergence failure — replicate remains NaN
                            % and is excluded from CI percentiles below
                        end
 
                        % Progress indicator: one dot per 100 iterations
                        if mod(b, 100) == 0, fprintf('.'); end
 
                    end
                    fprintf(' Done.\n');
 
                    % 95% percentile CI (NaN replicates excluded)
                    validBoots = bootICCs(~isnan(bootICCs));
                    ci_low     = prctile(validBoots, 2.5);
                    ci_high    = prctile(validBoots, 97.5);
 
                    %  STEP 5 — APPEND RESULT ROW
                    resultRow = table( ...
                        bl, {currentRef}, {currentAlgo}, {targetMetric}, ...
                        var_between, var_within, ...
                        (var_between / total_var) * 100, ...
                        (var_within  / total_var) * 100, ...
                        icc_single, ci_low, ci_high, ...
                        dependability_g, mean(countsPerPerson), ...
                        'VariableNames', { ...
                            'Baseline', 'Reference', 'Algorithm', 'Measurement', ...
                            'Var_Between', 'Var_Within', ...
                            'Var_Between_Pct', 'Var_Within_Pct', ...
                            'ICC_SingleTrial', 'CI_95_Low', 'CI_95_High', ...
                            'Dependability_G', 'Avg_Trials'});
 
                    allResults = [allResults; resultRow]; %#ok<AGROW>
 
                catch ME
                    fprintf('    !! Error in %s: %s\n', targetMetric, ME.message);
                end
 
            end % metrics loop
        end % reference loop
    end % baseline loop
end % algorithm loop
 
%  SECTION 4 — SAVE RESULTS
% ============================
 
if isempty(allResults)
    warning('No results were collected. Check file paths and condition filter.');
else
    writetable(allResults, outFile, 'Delimiter', ';');
 
    fprintf('Reliability analysis complete.\n');
    fprintf('  CSV results : %s\n', outFile);
