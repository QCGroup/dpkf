% clc
% close
% clear

%% Initialization

% Number of trials
numTrials = 200;

% Hovercraft dynamics
A = [cos(pi/18) -sin(pi/18); sin(pi/18) cos(pi/18)];
C = [1 1];
N = 200;
P0 = 1*eye(2);
mu0 = [1 1]';

% Process noise (zero-mean Gaussian)
Q = 5e-2*eye(2);
sqrtQ = sqrtm(Q);

% Measurement noise (exponential)
lambda = sqrt(1/3);

% Proposed estimator initialization
mu_proposed(:, 1) = mu0;
P_proposed{1} = P0;

% Kalman filter initialization
mu_kf(:, 1) = mu0;
P_kf{1} = P0;
R = 1/lambda^2;
muV = 1/lambda;

% Particle filter initialization
numParticles = 1000;
mu_pf(:, 1) = mu0;
P_pf{1} = P0;

% MCKF initialization
kernel = 4;
epsilon = 1e-6;
mu_mckf(:, 1) = mu0;
P_mckf{1} = P0;

% Masreliez filter initialization
mu_mas(:, 1) = mu0;
P_mas{1} = P0;

%% Main loop

for trialIter = 1:numTrials
    
    disp(['iteration ', num2str(trialIter), ' out of ', num2str(numTrials)])
    
    % Initial condition
    rng(trialIter, "twister");
    mu_proposed(:, 1) = mu0;
    mu_kf(:, 1) = mu0;
    mu_pf(:, 1) = mu0;
    mu_mas(:, 1) = mu0;
    mu_gsf(:, 1) = mu0;
    X0 = mu0 + sqrtm(P0)*randn(2, 1);
    X(:, 1) = X0;
    pf = particleFilter(@stateFcn, @LikelihoodFcs);
    initialize(pf, numParticles, mu0, P0);  

    % Measurement noise
    v = exprnd(muV, N+1, 1);

    % Process noise
    w = sqrtQ*randn(2, N+1);
    
    % Simulate the system and estimators
    for t = 1:N
    
        % Simulate the system
        X(:, t+1) = A*X(:, t) + w(:, t);
    
        % Recieve a measurement
        y(:, t+1) = C*X(:, t+1) + v(t+1);
        
        % Soluting using Kalman. Assuming both w and v are Gaussian
        tic;
        [mu_kf(:, t+1), P_kf{t+1}] = KF(mu_kf(:, t), P_kf{t}, muV, y(:, t+1), Q, R, A, C);
        time_kf{trialIter, 1}(t) = toc;
    
        % Soluting using the proposed estimator
        tic;
        [mu_proposed(:, t+1), P_proposed{t+1}] = poposedEstimator(mu_proposed(:, t), P_proposed{t}, lambda, C, A, y(:, t+1), Q);
        time_prposed{trialIter, 1}(t) = toc;
                                                                   
        % Solution using particle filter
        tic;
        predict(pf, numParticles, A, sqrtQ);
        correct(pf, y(:, t+1), C, lambda);
        [mu_pf(:, t+1), P_pf{t+1}] = getStateEstimate(pf);
         time_pf{trialIter, 1}(t) = toc;
        
        % Solution using MCKF
        tic;
        [mu_mckf(:, t+1), P_mckf{t+1}, ~] = MCKF(A, mu_mckf(:, t), P_mckf{t}, C, y(:, t+1), Q, R, muV, kernel, epsilon);
        time_mckf{trialIter, 1}(t) = toc;

%         % Solution using Masrelirz filter
%         tic;
%         [mu_mas(:, t+1), P_mas{t+1}] = Masreliez(A, C, Q, y(:, t+1), lambda, mu_mas(:, t), P_mas{t}, 0);
%         time_mas{trialIter, 1}(t) = toc;
        
    end

    % Compute position RMSE
    for i=1:length(mu_proposed)
        se_proposed = 0;
        se_kf = 0;
        se_mckf = 0;
        se_pf = 0;
%         se_mas = 0;
        for j=1:i
            pos_proposed = sqrt(norm(mu_proposed(1:2, j)));
            pos_kf = sqrt(norm(mu_kf(1:2, j)));
            pos_mckf = sqrt(norm(mu_mckf(1:2, j)));
            pos_pf = sqrt(norm(mu_pf(1:2, j)));
%             pos_mas = sqrt(norm(mu_mas(1:2, j)));
            pos_true = sqrt(norm(X(1:2, j)));
            se_proposed = se_proposed + (pos_proposed - pos_true)^2;
            se_kf = se_kf + (pos_kf - pos_true)^2;
            se_mckf = se_mckf + (pos_mckf - pos_true)^2;
            se_pf = se_pf + (pos_pf - pos_true)^2;
%             se_mas = se_mas + (pos_mas - pos_true)^2;
        end
        rmse_proposed{trialIter, 1}(i) = sqrt(se_proposed/i);
        rmse_kf{trialIter, 1}(i) = sqrt(se_kf/i);
        rmse_mckf{trialIter, 1}(i) = sqrt(se_mckf/i);
        rmse_pf{trialIter, 1}(i) = sqrt(se_pf/i);
%         rmse_mas{trialIter, 1}(i) = sqrt(se_mas/i);
    end
end

%% Plots

for i=1:N+1
    tmp = cell2mat(rmse_proposed);
    rmse_proposed_mean(i) = mean(tmp(:, i));
    rmse_proposed_std(i) = std(tmp(:, i));

    tmp = cell2mat(rmse_kf);
    rmse_kf_mean(i) = mean(tmp(:, i));
    rmse_kf_std(i) = std(tmp(:, i));

    tmp = cell2mat(rmse_mckf);
    rmse_mckf_mean(i) = mean(tmp(:, i));
    rmse_mckf_std(i) = std(tmp(:, i));

    tmp = cell2mat(rmse_pf);
    rmse_pf_mean(i) = mean(tmp(:, i));
    rmse_pf_std(i) = std(tmp(:, i));
 
%     tmp = cell2mat(rmse_mas);
%     rmse_mas_mean(i) = mean(tmp(:, i));
%     rmse_mas_std(i) = std(tmp(:, i));
end

% z = 1;
% faceAlpha = 0.3;
% f = figure;
% f.Units = "points";
% f.Position = [1 1 245.7 140];
% hold on;
% box on;
% grid on;
% set(0, 'DefaultTextInterpreter', 'latex');
% set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
% set(0, 'DefaultLegendInterpreter', 'latex');
% set(0, 'DefaultColorbarTickLabelInterpreter', 'latex');
% set(0, 'DefaultAxesFontName', 'TimesNewRoman');
% set(0, 'DefaultTextFontName', 'TimesNewRoman');
% set(0, 'DefaultAxesFontSize', 8);
% set(0, 'DefaultTextFontSize', 10);
% set(0, 'DefaultLineLineWidth', 0.7);
% set(0, 'DefaultLineLineStyle', '-');
% plot(0:N, rmse_proposed_mean, 'Color', [0 0.4470 0.7410]);
% fill([0:N N:-1:0], [rmse_proposed_mean+z*rmse_proposed_std rmse_proposed_mean(end:-1:1)-z*rmse_proposed_std(end:-1:1)], [0 0.4470 0.7410],...
%     'FaceAlpha', faceAlpha, 'HandleVisibility', 'off', 'EdgeColor', 'none');
% plot(0:N, rmse_kf_mean, 'Color', [0.8500 0.3250 0.0980]);
% fill([0:N N:-1:0], [rmse_kf_mean+z*rmse_kf_std rmse_kf_mean(end:-1:1)-z*rmse_kf_std(end:-1:1)], [0.8500 0.3250 0.0980], ...
%     'FaceAlpha', faceAlpha, 'HandleVisibility', 'off', 'EdgeColor', 'none');
% plot(0:N, rmse_mckf_mean, 'k');
% fill([0:N N:-1:0], [rmse_mckf_mean+z*rmse_mckf_std rmse_mckf_mean(end:-1:1)-z*rmse_mckf_std(end:-1:1)], 'k', ...
%     'FaceAlpha', faceAlpha, 'HandleVisibility', 'off', 'EdgeColor', 'none');
% % plot(0:N, rmse_mas_mean, 'Color', [0.6350 0.0780 0.1840]);
% % fill([0:N N:-1:0], [rmse_mas_mean+z*rmse_mas_std rmse_mas_mean(end:-1:1)-z*rmse_mas_std(end:-1:1)], [0.6350 0.0780 0.1840], ...
% %     'FaceAlpha', faceAlpha, 'HandleVisibility', 'off', 'EdgeColor', 'none');
% plot(0:N, rmse_pf_mean, 'Color', [0.4660 0.6740 0.1880]);
% fill([0:N N:-1:0], [rmse_pf_mean+z*rmse_pf_std rmse_pf_mean(end:-1:1)-z*rmse_pf_std(end:-1:1)], [0.4660 0.6740 0.1880],...
%     'FaceAlpha', faceAlpha, 'HandleVisibility', 'off', 'EdgeColor', 'none');
% ylabel('Position RMSE');
% xlabel('Time step');
% lgd = legend('Proposed', 'KF', 'MCKF', 'Masreliez', 'PF');
% lgd.Box = "off";
% lgd.ItemTokenSize = [4, 12, 12, 12];
% lgd.Location = "northoutside";
% lgd.Orientation = "horizontal";
% exportgraphics(f, 'fig_exponential.pdf', 'ContentType','vector');


%% Show the mean and std of rmse

% clc
% disp('Simulation for Gaussian process noise and Gaussian mixture measurement noise.');
% disp(['RMSE of the KF has mean ', num2str(rmse_kf_mean(end)), ' and std ', num2str(rmse_kf_std(end)), '.']);
% disp(['RMSE of the MCKF has mean ', num2str(rmse_mckf_mean(end)), ' and std ', num2str(rmse_mckf_std(end)), '.']);
% % disp(['RMSE of the Masreliez has mean ', num2str(rmse_mas_mean(end)), ' and std ', num2str(rmse_mas_std(end)), '.']);
% disp(['RMSE of the Proposed method has mean ', num2str(rmse_proposed_mean(end)), ' and std ', num2str(rmse_proposed_std(end)), '.']);
% disp(['RMSE of the PF has mean ', num2str(rmse_pf_mean(end)), ' and std ', num2str(rmse_pf_std(end)), '.'])

%% Computation time
for i=1:N
    tmp = cell2mat(time_prposed);
    time_proposed_mean(i) = mean(tmp(:, i));
    time_proposed_std(i) = std(tmp(:, i));

%     tmp = cell2mat(time_mas);
%     time_mas_mean(i) = mean(tmp(:, i));
%     time_mas_std(i) = std(tmp(:, i));

    tmp = cell2mat(time_kf);
    time_kf_mean(i) = mean(tmp(:, i));
    time_kf_std(i) = std(tmp(:, i));

    tmp = cell2mat(time_pf);
    time_pf_mean(i) = mean(tmp(:, i));
    time_pf_std(i) = std(tmp(:, i));

    tmp = cell2mat(time_mckf);
    time_mckf_mean(i) = mean(tmp(:, i));
    time_mckf_std(i) = std(tmp(:, i));

end

% disp('Simulation for Gaussian process noise and bimodal measurement noise.');
% disp('Computation time of KF is 1');
% disp(['Computation time of MCKF is ', num2str(mean(time_mckf_mean) / mean(time_kf_mean))]);
% % disp(['Computation time of Masreliez is ', num2str(mean(time_mas_mean) / mean(time_kf_mean))]);
% disp(['Computation time of Proposed is ', num2str(mean(time_proposed_mean) / mean(time_kf_mean))]);
% disp(['Computation time of PF is ', num2str(mean(time_pf_mean) / mean(time_kf_mean))]);
% fprintf('\n')

%% Save the results
clear f
if ~exist('results', 'dir')
    mkdir('results');
end
save(fullfile('results', 'workspace.mat'));

%% Functions used for particel filter

% State transition function
function out = stateFcn(particles, numParticles, A, sqrtQ)
    for i=1:numParticles
        particles(:, i) = A*particles(:, i);
    end
    w = sqrtQ*randn(2, numParticles);
    out = particles + w;
end

% Measurement function
function out = LikelihoodFcs(predictedParticels, y_t, C, lambda)
    residuals = y_t - C*predictedParticels;
    L = zeros(size(residuals));
    ok = residuals >= 0;
    L(ok) = lambda .* exp(-lambda .* residuals(ok));
    out = L;
end