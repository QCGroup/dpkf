function exportTableII(baseDir)
% EXPORTTABLEII_MATLAB (Option A - Screenshot)
%   Creates Table II (RMSE ± s.d. and normalized Time) from saved results
%   and exports a screenshot as TableII.pdf (bitmap, no LaTeX).

if nargin < 1, baseDir = pwd; end

noise_types  = {'skewedNormal','bimodal','gamma','impulsive', ...
                'cauchy','betaPrime','exponential','levy'};
noise_labels = {'Skewed normal (a)','Bimodal (b)','Gamma (c)','Impulsive (d)', ...
                'Cauchy (e)','Beta prime (f)','Exponential (g)','Lévy (h)'};
estimators   = {'KF','MCKF','Masreliez','Proposed','PF'};

E = numel(estimators);
M = numel(noise_types);

% Preallocate string arrays
RMSE = strings(E,M);
TIME = strings(E,M);

for j = 1:M
    ws = fullfile(baseDir, noise_types{j}, 'results', 'workspace.mat');
    if ~isfile(ws), continue, end
    S = load(ws);

    % Check availability
    hasKF   = isfield(S,'rmse_kf_mean');
    baselineVar = 'time_proposed';
    if isfield(S,'time_prposed'), baselineVar = 'time_prposed'; end
    if hasKF, baselineVar = 'time_kf'; end

    % --- RMSE ---
    RMSE(1,j) = valOrNA(S,'rmse_kf_mean','rmse_kf_std');
    RMSE(2,j) = valOrNA(S,'rmse_mckf_mean','rmse_mckf_std');
    RMSE(3,j) = valOrDiverged(S,'rmse_mas_mean','rmse_mas_std');
    RMSE(4,j) = valOrNA(S,'rmse_proposed_mean','rmse_proposed_std');
    RMSE(5,j) = valOrNA(S,'rmse_pf_mean','rmse_pf_std');

    % --- Time ---
    if hasKF
        TIME(1,j) = "1.00";
    else
        TIME(1,j) = "—";
    end
    TIME(2,j) = timeRatio(S,'time_mckf',baselineVar);
    TIME(3,j) = timeRatio(S,'time_mas', baselineVar);
    TIME(4,j) = timeRatio(S,pickPropVar(S),baselineVar,true);
    TIME(5,j) = timeRatio(S,'time_pf', baselineVar);
end

% ==== Build two block tables ====
blocks = {1:4, 5:8};   % first 4 noises, last 4 noises
pngPath = fullfile('figs','TableII.png');

f = figure('Visible','on','Units','normalized','Position',[0.1 0.1 0.8 0.8]);

for b = 1:2
    cols = blocks{b};
    dataBlock = cell(E,1+2*numel(cols));
    dataBlock(:,1) = estimators';
    colnames = {'Estimator'};
    for jj = 1:numel(cols)
        j = cols(jj);
        colnames{end+1} = [noise_labels{j} ' RMSE ± s.d.'];
        colnames{end+1} = [noise_labels{j} ' Time'];
        for e = 1:E
            dataBlock{e,2*jj}   = char(RMSE(e,j));  % string → char
            dataBlock{e,2*jj+1} = char(TIME(e,j));
        end
    end
    % Ensure char
    dataBlock = cellfun(@char, dataBlock, 'UniformOutput', false);

    uitable('Parent',f, ...
        'Data',dataBlock, ...
        'ColumnName',colnames, ...
        'RowName',[], ...
        'Units','normalized', ...
        'Position',[0, 1-b*0.5, 1, 0.5]); % top half then bottom half
end

drawnow; pause(0.5);  % let MATLAB render

% Capture screenshot
frame = getframe(f);
im = frame.cdata;

% Save PNG and PDF (bitmap)
imwrite(im, pngPath);
% imwrite(im, pdfPath);

close(f)
% fprintf('Saved %s and %s\n', pngPath, pdfPath);

end

% ===== Helper functions =====
function s = valOrNA(S,meanName,stdName)
    if isfield(S,meanName) && isfield(S,stdName) ...
            && isfinite(S.(meanName)(end)) && isfinite(S.(stdName)(end))
        s = sprintf('%.3f ± %.3f',S.(meanName)(end),S.(stdName)(end));
    else
        s = "N/A";
    end
end

function s = valOrDiverged(S,meanName,stdName)
    if isfield(S,meanName) && isfield(S,stdName)
        mu = S.(meanName)(end); sd = S.(stdName)(end);
        if isfinite(mu) && isfinite(sd)
            s = sprintf('%.3f ± %.3f',mu,sd);
            return
        end
    end
    s = "diverged";
end

function nm = pickPropVar(S)
    if isfield(S,'time_prposed'), nm='time_prposed'; else, nm='time_proposed'; end
end

function s = timeRatio(S,varName,baseVar,isProp)
    if nargin<4, isProp=false; end
    if ~isfield(S,varName) || ~isfield(S,baseVar)
        s="—"; return
    end
    N = size(cell2mat(S.(varName)),2);
    gm1 = gmMean(S.(varName),N);
    gm2 = gmMean(S.(baseVar),N);
    if ~isfinite(gm1) || ~isfinite(gm2) || gm2==0
        s="—";
    else
        if isProp && strcmp(varName,baseVar)
            s="1.00";
        else
            s=sprintf('%.2f',gm1/gm2);
        end
    end
end

function g = gmMean(cellTimes,N)
    M = cell2mat(cellTimes);
    v=zeros(1,N);
    for i=1:N
        c = M(:,i); c=c(c>0);
        if isempty(c), v(i)=NaN; else, v(i)=geo_mean(c); end
    end
    v=v(isfinite(v));
    if isempty(v), g=NaN; else, g=geo_mean(v); end
end
