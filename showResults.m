% This code uses the stored data to show the results and export their
% tables as PDFs.



%% main code

clc
close 
clear
baseDir = pwd;
numNoise = 8;
clc
for iter=1:numNoise
    clearvars -except iter baseDir numNoise
    noise_types = {'skewedNormal', 'bimodal', 'gamma', 'impulsive', 'cauchy', 'betaPrime', 'exponential', 'levy'};
    inner = noise_types{iter};
    resultsFile = fullfile(baseDir, inner, 'results/workspace.mat');
    load(resultsFile);
    clear noise_types
    noise_types = {'skewedNormal', 'bimodal', 'gamma', 'impulsive', 'cauchy', 'betaPrime', 'exponential', 'levy'};
    disp(['RMSE (steady state) results for ', noise_types{iter}, ' measurement noise.']);
    if ~exist('rmse_mas_mean', 'var') && ~exist('rmse_kf_mean', 'var')
        disp(['The Proposed method has mean ', num2str(rmse_proposed_mean(end)), ' and std ', num2str(rmse_proposed_std(end)), '.']);
        disp(['The PF has mean ', num2str(rmse_pf_mean(end)), ' and std ', num2str(rmse_pf_std(end)), '.']);
    elseif ~exist('rmse_mas_mean', 'var')
        disp(['The KF has mean ', num2str(rmse_kf_mean(end)), ' and std ', num2str(rmse_kf_std(end)), '.']);
        disp(['The MCKF has mean ', num2str(rmse_mckf_mean(end)), ' and std ', num2str(rmse_mckf_std(end)), '.']);
        disp(['The Proposed method has mean ', num2str(rmse_proposed_mean(end)), ' and std ', num2str(rmse_proposed_std(end)), '.']);
        disp(['The PF has mean ', num2str(rmse_pf_mean(end)), ' and std ', num2str(rmse_pf_std(end)), '.']);
    elseif ~exist('rmse_kf_mean', 'var')
        disp(['The Masreliez has mean ', num2str(rmse_mas_mean(end)), ' and std ', num2str(rmse_mas_std(end)), '.']);
        disp(['The Proposed method has mean ', num2str(rmse_proposed_mean(end)), ' and std ', num2str(rmse_proposed_std(end)), '.']);
        disp(['The PF has mean ', num2str(rmse_pf_mean(end)), ' and std ', num2str(rmse_pf_std(end)), '.']);
    else
        disp(['The KF has mean ', num2str(rmse_kf_mean(end)), ' and std ', num2str(rmse_kf_std(end)), '.']);
        disp(['The MCKF has mean ', num2str(rmse_mckf_mean(end)), ' and std ', num2str(rmse_mckf_std(end)), '.']);
        disp(['The Masreliez has mean ', num2str(rmse_mas_mean(end)), ' and std ', num2str(rmse_mas_std(end)), '.']);
        disp(['The Proposed method has mean ', num2str(rmse_proposed_mean(end)), ' and std ', num2str(rmse_proposed_std(end)), '.']);
        disp(['The PF has mean ', num2str(rmse_pf_mean(end)), ' and std ', num2str(rmse_pf_std(end)), '.']);
    end
    disp(['Computation time (average) results for ', noise_types{iter}, ' measurement noise.']);
    if ~exist('rmse_mas_mean', 'var') && ~exist('rmse_kf_mean', 'var')
        for i=1:N
            tmp = cell2mat(time_prposed);
            time_proposed_mean(i) = geo_mean(tmp(:, i));
            tmp = cell2mat(time_pf);
            time_pf_mean(i) = geo_mean(tmp(:, i));
        end

        disp(['Computation time of Proposed is ', num2str(geo_mean(time_proposed_mean) / geo_mean(time_proposed_mean))]);
        disp(['Computation time of PF is ', num2str(geo_mean(time_pf_mean) / geo_mean(time_proposed_mean))]);
    elseif ~exist('rmse_mas_mean', 'var')
        for i=1:N
            tmp = cell2mat(time_prposed);
            time_proposed_mean(i) = geo_mean(tmp(:, i));
            tmp = cell2mat(time_kf);
            time_kf_mean(i) = geo_mean(tmp(:, i));
            tmp = cell2mat(time_pf);
            time_pf_mean(i) = geo_mean(tmp(:, i));
            tmp = cell2mat(time_mckf);
            time_mckf_mean(i) = geo_mean(tmp(:, i));
        end
        disp('Computation time of KF is 1');
        disp(['Computation time of MCKF is ', num2str(geo_mean(time_mckf_mean) / geo_mean(time_kf_mean))]);
        disp(['Computation time of Proposed is ', num2str(geo_mean(time_proposed_mean) / geo_mean(time_kf_mean))]);
        disp(['Computation time of PF is ', num2str(geo_mean(time_pf_mean) / geo_mean(time_kf_mean))]);
    elseif ~exist('rmse_kf_mean', 'var')
        for i=1:N
            tmp = cell2mat(time_prposed);
            time_proposed_mean(i) = geo_mean(tmp(:, i));
            tmp = cell2mat(time_mas);
            time_mas_mean(i) = geo_mean(tmp(:, i));
            tmp = cell2mat(time_pf);
            time_pf_mean(i) = geo_mean(tmp(:, i));
        end
        disp(['Computation time of Masreliez is ', num2str(geo_mean(time_mas_mean) / geo_mean(time_proposed_mean))]);
        disp(['Computation time of Proposed is ', num2str(geo_mean(time_proposed_mean) / geo_mean(time_proposed_mean))]);
        disp(['Computation time of PF is ', num2str(geo_mean(time_pf_mean) / geo_mean(time_proposed_mean))]);
    else
        for i=1:N
            tmp = cell2mat(time_prposed);
            time_proposed_mean(i) = geo_mean(tmp(:, i));
            tmp = cell2mat(time_mas);
            time_mas_mean(i) = geo_mean(tmp(:, i));
            tmp = cell2mat(time_kf);
            time_kf_mean(i) = geo_mean(tmp(:, i));
            tmp = cell2mat(time_pf);
            time_pf_mean(i) = geo_mean(tmp(:, i));
            tmp = cell2mat(time_mckf);
            time_mckf_mean(i) = geo_mean(tmp(:, i));
        end
        disp('Computation time of KF is 1');
        disp(['Computation time of MCKF is ', num2str(geo_mean(time_mckf_mean) / geo_mean(time_kf_mean))]);
        disp(['Computation time of Masreliez is ', num2str(geo_mean(time_mas_mean) / geo_mean(time_kf_mean))]);
        disp(['Computation time of Proposed is ', num2str(geo_mean(time_proposed_mean) / geo_mean(time_kf_mean))]);
        disp(['Computation time of PF is ', num2str(geo_mean(time_pf_mean) / geo_mean(time_kf_mean))]);
    end
    fprintf('\n\n\n');
end
clear