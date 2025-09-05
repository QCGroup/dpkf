clc
clear
close

n = 200;

%% beta prime
x_bp = linspace(0,6,n);
alpha = 2;
beta = fzero(@(b) alpha*(alpha + b - 1) - 3*(b - 2)*(b - 1)^2, 3);
p_bp = betaprime_pdf(x_bp, alpha, beta);

%% bimodal measurement noise
x_bm = linspace(-5,5,n);
p_bm = 0.4 * normpdf(x_bm,-1.8,sqrt(0.9)) + 0.6 * normpdf(x_bm,1.2,sqrt(0.8));

%% Cauchy
x_c = linspace(-5,5,n);
gamma = 1;
p_c = cauchy_pdf(x_c, gamma);

%% Exponential
x_e = linspace(0,5,n);
lambda = sqrt(1/3);
p_e = exponential_pdf(x_e, lambda);

%% Gamma
x_g = linspace(0,8,n);
alpha = 2;
theta = sqrt(3/2);
p_g = gamma_pdf(x_g, alpha, theta);

%% Impulsive noise
x_im = linspace(-20,20,10*n);
p_im = 0.1 * normpdf(x_im,0,sqrt(25)) + 0.9 * normpdf(x_im,0,sqrt(0.5556));

%% Skewed normal measurement noise
x_sn = linspace(-5,10,n);
p_sn = SN(x_sn,-2.006269,2.6505,3);

%% Levy
x_l = linspace(1,10,n);
mu = 1;
c = 3;
p_l = levy_pdf(x_l, mu, c);

%% Tiled plot 
close
figure('Name','1', 'Units', 'points', 'Position', [1 1 505.89 120]);
clf
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultColorbarTickLabelInterpreter', 'latex');

t = tiledlayout(2, 4, "TileSpacing", "loose", "Padding", "none");
fontSize = 7;
lineWidth = 1;

nexttile
plot(x_sn, p_sn, 'Color', "#4575b4", 'LineWidth', lineWidth, 'LineStyle', '-');
grid on;
ylim('tight');
xlim('tight');
% xticks(min(x_sn):1:max(x_sn))
ylim([0 max(ylim)])
set(gca, 'FontName','Times New Roman')
toprighttext('(a)');

nexttile
plot(x_bm, p_bm, 'Color', "#4575b4", 'LineWidth', lineWidth, 'LineStyle', '-');
grid on;
ylim('tight');
xlim('tight');
yticks([0 0.1 0.2 0.3]);
ylim([0 max(ylim)])
set(gca, 'FontName','Times New Roman')
toprighttext('(b)')

nexttile
plot(x_g, p_g, 'Color', "#4575b4", 'LineWidth', lineWidth, 'LineStyle', '-');
grid on;
ylim("tight");
xlim('tight');
% xticks([-1, 0, 1]);
ylim([0 max(ylim)]);
set(gca, 'FontName','Times New Roman')  
toprighttext('(c)')

nexttile
plot(x_im, p_im, 'Color', "#4575b4", 'LineWidth', lineWidth, 'LineStyle', '-');
ax = gca;
grid on;
set(gca, 'YScale', 'log');  
ylim('tight');
xlim('tight');
yticks(10.^[-6, -4, -2, 0]);
ylim([10^(-6), 10^0]);
ax.YMinorGrid = 'off';
ax.YGrid = 'on';
ax.YMinorTick = 'off';
toprighttext('(d)')
set(gca, 'FontName','Times New Roman');

nexttile
plot(x_c, p_c, 'Color', "#4575b4", 'LineWidth', lineWidth, 'LineStyle', '-');
grid on;
ylim("tight");
xlim('tight');
% xticks([-1, 0, 1]);
ylim([0 max(ylim)]);
set(gca, 'FontName','Times New Roman')  
toprighttext('(e)')

nexttile
plot(x_bp, p_bp, 'Color', "#4575b4", 'LineWidth', lineWidth, 'LineStyle', '-');
grid on;
ylim("tight");
xlim('tight');
% xticks([-1, 0, 1]);
ylim([0 max(ylim)]);
set(gca, 'FontName','Times New Roman')  
toprighttext('(f)')

nexttile
plot(x_e, p_e, 'Color', "#4575b4", 'LineWidth', lineWidth, 'LineStyle', '-');
grid on;
ylim("tight");
xlim('tight');
% xticks([-1, 0, 1]);
ylim([0 max(ylim)]);
set(gca, 'FontName','Times New Roman')  
toprighttext('(g)')

nexttile
plot(x_l, p_l, 'Color', "#4575b4", 'LineWidth', lineWidth, 'LineStyle', '-');
grid on;
ylim("tight");
xlim('tight');
% xticks([-1, 0, 1]);
ylim([0 max(ylim)]);
set(gca, 'FontName','Times New Roman')  
toprighttext('(h)')



%% Export the figure
exportgraphics(t, 'fig_pdf_samples.pdf', 'ContentType','vector');

%% Helper functions 

function toprighttext(txt)
    % Get axis limits
    xlimVals = xlim;
    ylimVals = ylim;

    % Small margin (e.g., 5% of axis range)
    xMargin = 0.05 * range(xlimVals);
    yScale = get(gca, 'YScale');
    switch yScale
        case 'linear'
            yMargin = 0.07 * range(ylimVals);
        case 'log'
            yMargin = 10^(0.07 * range(log10(ylimVals)));
    end

    % Define position near the upper right corner
    xText = xlimVals(2) - xMargin;
    switch yScale
        case 'linear'
            yText = ylimVals(2) - yMargin;
        case 'log'
            yText = ylimVals(2)/yMargin;
    end

    % Place text
    text(xText, yText, txt, ...
         'HorizontalAlignment', 'right', ...
         'VerticalAlignment', 'top',...
         'BackgroundColor', 'white', ...
         'FontWeight', 'bold', ...
         'FontSize', 8, ...
         'Margin',2);
end



%% Distribtuion functions 

function p = betaprime_pdf(x, alpha, beta)
    p = zeros(size(x));
    
    % Valid support: x > 0
    valid = x > 0;
    if any(valid)
        coeff = 1 / beta_func(alpha, beta);
        p(valid) = coeff .* (x(valid).^(alpha - 1)) .* ((1 + x(valid)).^(-(alpha + beta)));
    end
end

function B = beta_func(a, b)
    B = gamma(a) * gamma(b) / gamma(a + b);
end

function p = SN(x,xi,s,alpha)
    p = 2/s.*normpdf((x-xi)/s).*normcdf(alpha*(x-xi)/s);
end

function p = cauchy_pdf(x, gamma)
    p = 1 ./ (pi * gamma * (1 + (x./gamma).^2));
end

function p = exponential_pdf(x, lambda)
    p = zeros(size(x));
    
    % Valid support: x >= 0
    valid = x >= 0;
    if any(valid)
        p(valid) = lambda .* exp(-lambda .* x(valid));
    end
end

function p = gamma_pdf(x, alpha, theta)
    p = zeros(size(x));
    
    % Valid support: x >= 0
    valid = x >= 0;
    if any(valid)
        coeff = 1 / (gamma(alpha) * theta^alpha);
        p(valid) = coeff .* (x(valid).^(alpha - 1)) .* exp(-x(valid) ./ theta);
    end
end

function p = levy_pdf(x, mu, c)    
    p = zeros(size(x));

    % Valid support: x > mu
    valid = x > mu;
    if any(valid)
        coeff = sqrt(c/(2*pi));
        p(valid) = coeff .* exp(-c ./ (2*(x(valid) - mu))) ./ ((x(valid) - mu).^(3/2));
    end
end
