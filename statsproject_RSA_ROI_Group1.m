%% SoSe23: Probabilistic and Statistical Modelling II
% Decoding Project, Group 1
% Student Name: Lok Yan Lam
% Last Update on: 21.08.2023
%% RSA based on cross-validated Euclidean distance

% Prerequisites for running this scripts:
% 
% % beta values from an SPM.mat file: check if beta_000?.nii files and
% % SPM.mat file exist and set up path for the folder that contains all 
% % experimental data at datadir
% % 
% % open a subfolder named 'results' in datadir
% %
% % SPM12: set up path at spmdir
% %
% % The Decoding Toolbox: set up path at tdtdir

% This script runs a model-based RSA in each ROI by correlating
% stimuli-based RDMs with a 6x6 neural RDM, constrcuted by leave-one-run-out
% sixfold cross-validated Euclidean distances between six different
% content-specific stimuli. Decoding accuracies are subsequently tested
% against chance level with one-sample t-tests.

%% Housekeeping

clear all;
clc;
%% Set paths

datadir = '/Volumes/Transcend/shared_data';
spmdir = '/Volumes/Transcend/MCNB - spm/spm12';
tdtdir = '/Volumes/Transcend/MCNB - spm/tdt_3.999F/decoding_toolbox';

addpath(datadir, genpath(spmdir), genpath(tdtdir));

%% Create and plot stimuli-based RDMs

% prepare stimuli names for plotting

row_labels = {'StimPress', 'StimFlutt', 'StimVibro', 'ImagPress', 'ImagFlutt', 'ImagVibro'};

% 1. Imagery vs. Perception RDM

diff_cat = ones(3,3); % high distances between different categories
same_cat = zeros(3,3); % low distances between same categories

% concatenate the matrices to create a 6x6 RDM
rdm_upper = cat(2,same_cat,diff_cat);
rdm_lower = cat(2,diff_cat,same_cat);
rdm_imag_per = cat(1,rdm_upper,rdm_lower); 

% Plot
rdm_imag_per_jpg = figure;
imagesc(rdm_imag_per);
title('RDM: Imagery vs. Perception');
yticklabels(row_labels);
xticklabels(row_labels);
colorbar;

% Print distance values on each cell in the RDM
for i = 1:size(rdm_imag_per, 1)
    for j = 1:size(rdm_imag_per, 2)
        text(j, i, num2str(rdm_imag_per(i, j), '%.2f'), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end

axis equal tight;

% save RDM
saveas(rdm_imag_per_jpg,fullfile(datadir, 'results', 'rdm_imag_per.jpg'));



% 2. Press vs. Flutter vs. Vibro RDM

V1 = ones(1,6); % vector of 6 ones 
M1 = diag(V1) + diag(V1(1:3), 3) + diag(V1(1:3),-3); % create a 6x6 matrix with a diagonal of ones and off-diagonal elements
% that are also ones, but shifted by three positions both to the right and left.
rdm = ~M1; % flip the ones and zeros in the matrix
rdm_contspec = double(rdm); % create RDM by converting booleans to numerical values

% Plot
rdm_contspec_jpg = figure;
imagesc(rdm_contspec);
title('RDM: Press vs. Flutter vs. Vibro');
yticklabels(row_labels);
xticklabels(row_labels);
colorbar;

% Print distance values on each cell in the RDM
for i = 1:size(rdm_contspec, 1)
    for j = 1:size(rdm_contspec, 2)
        text(j, i, num2str(rdm_contspec(i, j), '%.2f'), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end

axis equal tight;

% save RDM
saveas(rdm_contspec_jpg,fullfile(datadir, 'results', 'rdm_contspec.jpg'));




% 3. 6-condition RDM

V2 = ones(1,3); % a vector of 3 ones
M2 = diag(V2); % create a 3x3 matrix with a diagonal of ones
rdm = ~M2; % flip the ones and zeros in the matrix
stim_rdm = double(rdm); % create within-category RDM by converting booleans to numerical values

imag_rdm = stim_rdm + 0.5; % introduce distance = 0.5 between the categories imagery and perception

% concatenate the matrices to create a 6x6 RDM
rdm_upper = cat(2,stim_rdm, imag_rdm);
rdm_lower = cat(2,imag_rdm, stim_rdm);
rdm_sixcond = cat(1,rdm_upper, rdm_lower);

% Plot
rdm_sixcond_jpg = figure;
imagesc(rdm_sixcond);
title('RDM: All six conditions');
yticklabels(row_labels);
xticklabels(row_labels);
colorbar;

% Print distance values on each cell in the RDM
for i = 1:size(rdm_sixcond, 1)
    for j = 1:size(rdm_sixcond, 2)
        text(j, i, num2str(rdm_sixcond(i, j), '%.2f'), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end

axis equal tight;

% save RDM
saveas(rdm_sixcond_jpg, fullfile(datadir, 'results', 'rdm_sixcond.jpg'));

%% Create neural RDMs based on leave-one-run-out sixfold cross-validated Euclidean distances & perform subject-level analysis

% list all ROIs to be investigated
rois = {
    '2CONJ_left_IPL_uncorr',...
    '3CONJ_SMA_uncorr',...
    'rPSC_2_TR50_right_CUT_Stim_vs_Null',...
    'rPSC_3b_TR50_right_CUT_Stim_vs_Null',...
    'rPSC_1_TR50_right_CUT_Stim_vs_Null'
    };

accuracies = []; % init empty array for storing decoding accuracies

for subj = 1:10

    % Set up configurations for analysis
    cfg = decoding_defaults;
    cfg.analysis = 'ROI'; % for ROI analysis

    % Create output directory for subject if it does not already exist
    subj_name = sprintf('sub-%03d', subj);
    OutfolderName = fullfile(datadir, 'results', subj_name);

    if ~exist(OutfolderName, 'dir')
        mkdir(OutfolderName);
    end

    for r = 1:length(rois)

        % Create output directory for each ROI under each subject if it does not already exist
        RoifolderName = fullfile(datadir, 'results', subj_name, rois{r});

        if ~exist(RoifolderName, 'dir')
            mkdir(RoifolderName);
        end

        % output directory
        cfg.results.dir = RoifolderName;

        % Set the filepath where SPM.mat and all related betas are
        beta_loc = fullfile(datadir, subj_name, '1st_level_good_bad_Imag');

        % Set the filename of brain mask
        cfg.files.mask = [datadir '/rois/' rois{1,r} '.nii'];

        % Set the label names to the regressor names for analysis
        labelnames = {'StimPress', 'StimFlutt', 'StimVibro', 'ImagPress', 'ImagFlutt', 'ImagVibro'};

        % each unique regressor name will receive a unique label (1 to n for n classes)
        labels = [1 2 3 4 5 6];

        % set everything to calculate (dis)similarity estimates
        cfg.decoding.software = 'distance'; % the difference to 'similarity' is that this averages across data with the same label
        cfg.decoding.method = 'classification'; % this is more a placeholder
        cfg.decoding.train.classification.model_parameters = 'cveuclidean'; % cross-validated Euclidean after noise normalization

        % Averages across cross-validation iterations
        cfg.results.output = 'other_average';

        % Extract all beta names and corresponding run numbers from the SPM.mat
        regressor_names = design_from_spm(beta_loc);

        % Extract all information for the cfg.files structure
        cfg = decoding_describe_data(cfg,labelnames,labels,regressor_names,beta_loc);

        % Create a design in which data is used to perform RSA with cross-validation
        cfg.design = make_design_similarity_cv(cfg);

        % Run decoding
        results = decoding(cfg);

        % load cross-validated neural RDM
        neural_rdm_file = load(fullfile(datadir, 'results', subj_name, rois{r}, 'res_other_average.mat'));
        cv_neural_rdm = neural_rdm_file.results.other_average.output{1,1};

        % Subject-level analysis: Compute decoding accuracies for each stimuli-based RDM
        accuracies(1,r,subj) = corr(rdm_imag_per(:), cv_neural_rdm(:), 'Type', 'Spearman');
        accuracies(2,r,subj) = corr(rdm_contspec(:), cv_neural_rdm(:), 'Type', 'Spearman');
        accuracies(3,r,subj) = corr(rdm_sixcond(:), cv_neural_rdm(:), 'Type', 'Spearman');

    end

end

% save decoding accurarcies
save(fullfile(datadir, 'results', 'accuracies.mat'), 'accuracies');


%% Group-level analysis

% load accuracy folder
data = load(fullfile(datadir, 'results', 'accuracies.mat'));
accuracies = data.accuracies;

% Define subject number (for standard error later)
subj_num = 10;

% Init an empty array to store the results for each hypothesis
stats_results = cell(3, length(rois));

% Define the false discovery rate significance level
fdr_threshold = 0.05;

% Init arrays to store the FDR-corrected significance level and
% significance markers for plotting asterisks later
fdr_corr_threshold = [];
significant_markers = zeros(3, length(rois));

% Init array to store p-values for FDR
p_values = zeros(3, length(rois));

% One-sample t-tests for each hypothesis in each ROI

for h = 1:3

    for r = 1:length(rois)
        % Perform t-tests for each hypothesis h for the current ROI (r)
        [~, p, ~, stats] = ttest(accuracies(h, r, :), 0);

        % Store the results in the cell array
        stats_results{h, r} = [stats.tstat, p];

        % Store the p-values for correction
        p_values(h,r) = p;

    end

    % Correct the significance level for multiple comparisons using the false discovery rate
    [is_significant, crit_p] = fdr_bh(p_values(h,:), fdr_threshold, 'pdep');

    % Store the FDR-corrected significance level and significance markers
    fdr_corr_threshold(h) = crit_p;
    significant_markers(h, :) = is_significant;
end

%% Plot results from group-level analysis

% Init the figure
plot = figure;

% Init arrays to store the mean_accuracies and std_errors in the required format for plotting bar chart

mean_accuracies = zeros(length(rois), 3);
std_errors = zeros(length(rois), 3);

% Prepare mean accurucies and standard error for plotting bar chart
for r = 1:length(rois)
    % Extract the accuracies and significant markers for the current ROI
    roi_accuracies = squeeze(accuracies(:, r, :)); % remove r (singleton) dimension
    roi_significant = significant_markers(:, r);
    
    % Calculate the mean and standard error of accuracies and store results
    % in the correponding arrays

    mean_accuracies(r,:) = (mean(roi_accuracies, 2))'; % mean accuracy for each hypothesis
    
    std_errors(r,:) = (std(roi_accuracies, 0, 2) / sqrt(subj_num))';
       
    % Create x-positions for the bars
    x_values = 1:length(rois);

end
    
    % Plot bar chart
    b = bar(mean_accuracies, 'grouped'); 
    b(2).FaceColor = [0.4660 0.6740 0.1880]; % change red bars to green for aesthetic reasons

    hold on;

    % Prepare errorbars

    [ngroups, nbars] = size(mean_accuracies); % get number of groups and bars

    % Get the x coordinates of the bars
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = b(i).XEndPoints;
    end

    % Plot the errorbars

    errorbar(x',mean_accuracies, std_errors,'k','linestyle','none');

    % Plot asterisks for significant results

    for h = 1:3
        for r = 1:length(rois)
            if significant_markers(h,r) % if the current hypothesis for the current ROI is significant
               asterisk_position = x(h,r); % extract the x-coordinate of the bar

               % Plot asterisk 0.05 higher than the upper limit of the
               % errorbar

               y = mean_accuracies(r,h) + std_errors(r,h);
               text(asterisk_position, y + 0.05, '*', 'FontSize', 16, 'HorizontalAlignment', 'center');
            end
        end
    end
    

% Customise the plot (labels, legend, etc.)
ylim([-0.2, 1.1])
xlabel('ROI');
ylabel('Mean Decoding Accuracy (Spearman Correlation)');
legend('Imag vs. Stim', 'Press vs. Flutter vs. Vibro', 'All six conditions');

% Set x-axis labels for ROIs
roi_labels = {'l IPL', 'SMA', 'r BA2', 'r BA 3b', 'r BA 1'}; 
xticks(1:length(rois));
xticklabels(roi_labels);

% Add title
title('Mean Decoding Accuracies for Different ROIs');

% Hold off to stop adding to the current figure
hold off;

% save the plot
saveas(plot, fullfile(datadir, 'results', 'group_level_stats.jpg'));