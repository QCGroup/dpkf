clc
clear
close
addpath('utils')

%% Plot Fig. 1 in the paper

plot_distributions;

%% Run the simulations or just show the results 

prompt = "If you want to show the results in Table II in the paper, press any key. If you want to run the entire simulations and store the results, enter 'y' . \n";
answer = input(prompt, "s");
if strcmp(answer, 'y') || strcmp(answer, 'Y') || strcmp(answer, 'yes') || strcmp(answer, 'Yes')
    disp('Running the entire simulations... Make sure to have the optimization toolbox installed on your machine...');
    pause(3);
    clc
    close
    clear
    baseDir = fullfile(pwd, 'utils');
    noise_types = {'betaPrime', 'bimodal', 'cauchy', 'exponential', 'gamma', 'impulsive', 'levy', 'skewedNormal'};
    numMeasNoise = length(noise_types);
    for iii=1:numMeasNoise
        pause(3);
        clc
        close
        clearvars -except iii numMeasNoise baseDir noise_types
        inner = noise_types{iii};
        mainFile = fullfile(baseDir, inner, 'main.m');
        if exist(mainFile, 'file')
            fprintf('Running: %s\n', mainFile);
            try
                run(mainFile);
            catch ME
                warning('Failed: %s\n  -> %s', mainFile, ME.message);
            end
        else
            warning('Missing file: %s', mainFile);
        end
    end
end
disp('showing the results from stored data.');
close 
clear
baseDir = fullfile(pwd, 'utils');
numNoise = 8;
clc
for iter=1:numNoise
    clearvars -except iter baseDir numNoise
    noise_types = {'skewedNormal', 'bimodal', 'gamma', 'impulsive', 'cauchy', 'betaPrime', 'exponential', 'levy'};
    inner = noise_types{iter};
    resultsFile = fullfile(baseDir, inner, 'results', 'workspace.mat');
    load(resultsFile);
    clear noise_types
    noise_types = {'skewedNormal', 'bimodal', 'gamma', 'impulsive', 'cauchy', 'betaPrime', 'exponential', 'levy'};
    disp(['RMSE results for ', noise_types{iter}, ' measurement noise.']);
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
exportTableII(baseDir)
disp("Table II has been generated and saved as .png file in the 'figs' folder");
clear

%% Plot Fig. 2 in the paper

baseDir = fullfile(pwd, 'utils');
resultsFile = fullfile(baseDir, 'impulsive', 'results', 'workspace1000.mat');
load(resultsFile);
z = 1.96;
faceAlpha = 0.3;
f = figure;
f.Units = "points";
f.Position = [1 1 245.7 140];
hold on;
box on;
grid on;
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultColorbarTickLabelInterpreter', 'latex');
set(0, 'DefaultAxesFontName', 'TimesNewRoman');
set(0, 'DefaultTextFontName', 'TimesNewRoman');
set(0, 'DefaultAxesFontSize', 8);
set(0, 'DefaultTextFontSize', 10);
set(0, 'DefaultLineLineWidth', 0.7);
set(0, 'DefaultLineLineStyle', '-');
plot(0:N, rmse_proposed_mean, 'Color', [0 0.4470 0.7410]);
fill([0:N N:-1:0], [rmse_proposed_mean+z*rmse_proposed_std rmse_proposed_mean(end:-1:1)-z*rmse_proposed_std(end:-1:1)], [0 0.4470 0.7410],...
    'FaceAlpha', faceAlpha, 'HandleVisibility', 'off', 'EdgeColor', 'none');
plot(0:N, rmse_kf_mean, 'Color', [0.8500 0.3250 0.0980]);
fill([0:N N:-1:0], [rmse_kf_mean+z*rmse_kf_std rmse_kf_mean(end:-1:1)-z*rmse_kf_std(end:-1:1)], [0.8500 0.3250 0.0980], ...
    'FaceAlpha', faceAlpha, 'HandleVisibility', 'off', 'EdgeColor', 'none');
plot(0:N, rmse_mckf_mean, 'k');
fill([0:N N:-1:0], [rmse_mckf_mean+z*rmse_mckf_std rmse_mckf_mean(end:-1:1)-z*rmse_mckf_std(end:-1:1)], 'k', ...
    'FaceAlpha', faceAlpha, 'HandleVisibility', 'off', 'EdgeColor', 'none');
plot(0:N, rmse_mas_mean, 'Color', [0.6350 0.0780 0.1840]);
fill([0:N N:-1:0], [rmse_mas_mean+z*rmse_mas_std rmse_mas_mean(end:-1:1)-z*rmse_mas_std(end:-1:1)], [0.6350 0.0780 0.1840], ...
    'FaceAlpha', faceAlpha, 'HandleVisibility', 'off', 'EdgeColor', 'none');
plot(0:N, rmse_pf_mean, 'Color', [0.4660 0.6740 0.1880], 'LineStyle', '-.');
fill([0:N N:-1:0], [rmse_pf_mean+z*rmse_pf_std rmse_pf_mean(end:-1:1)-z*rmse_pf_std(end:-1:1)], [0.4660 0.6740 0.1880],...
    'FaceAlpha', faceAlpha, 'HandleVisibility', 'off', 'EdgeColor', 'none');
ylabel('RMSE');
xlabel('Time step');
lgd = legend('Proposed', 'KF', 'MCKF', 'Masreliez', 'PF');
lgd.Box = "off";
lgd.ItemTokenSize = [8, 8, 8, 8];
lgd.Location = "northoutside";
lgd.Orientation = "horizontal";
xlim('tight');
ylim('tight');
exportgraphics(f, fullfile('figs', 'Fig2.pdf'), 'ContentType','vector');