function [g, G] = compute_gG(y_t, C, mu_tmp, P_tmp, gamma)
% Computes g(z) and G(z) using 2D quadrature
% z  : scalar measurement
% H  : 1x2 measurement matrix
% mu : 2x1 mean of Gaussian prior
% M  : 2x2 covariance of Gaussian prior

    % ----- Parameters for quadrature -----
     ns = 50;  % number of grid points per dimension (adjust as needed)
    L = 5;    % how many std devs to cover in each direction
    
    % Standard deviations from covariance
    stds = sqrt(diag(P_tmp));
    range1 = linspace(mu_tmp(1)-L*stds(1), mu_tmp(1)+L*stds(1), ns);
    range2 = linspace(mu_tmp(2)-L*stds(2), mu_tmp(2)+L*stds(2), ns);
    [X1, X2] = meshgrid(range1, range2);
    
    % Area element
    dx1 = range1(2)-range1(1);
    dx2 = range2(2)-range2(1);
    w = dx1*dx2;
    
    % Gaussian prior density at grid
    X = [X1(:), X2(:)]';
    prior = mvnpdf(X', mu_tmp', P_tmp);  % ns^2 x 1
    
    % Residuals u = z - Hx
    U = y_t - (C*X);  % 1 x ns^2
    
    % Evaluate noise pdf and derivatives (user fills these in)
    v0 = p(U, gamma);    % p_v(u)
    v1 = dp(U, gamma);   % derivative wrt u
    v2 = ddp(U, gamma);  % second derivative wrt u
    
    % Weighted sums (quadrature)
    P  = sum(v0(:) .* prior) * w;
    P1 = sum(v1(:) .* prior) * w;
    P2 = sum(v2(:) .* prior) * w;
    
    % Compute g and G
    g = - P1 / P;
    G = - (P2 / P - (P1 / P)^2);

    % --- Helper functions to be filled by user ---
    function y = p(u, gamma)
        % PDF of noise v at u
        y = 1 ./ (pi*gamma * (1 + (u./gamma).^2));
        % Fill this in
    end
    
    function y = dp(u, gamma)
        % First derivative of noise PDF wrt u
        y = - (2 .* u) ./ (pi * gamma^3 * (1 + (u./gamma).^2).^2);
        % Fill this in
    end
    
    function y = ddp(u, gamma)
        % Second derivative of noise PDF wrt u
        y = - (2 ./ (pi * gamma^3)) .* (1 - 3*(u./gamma).^2) ./ (1 + (u./gamma).^2).^3;
        % Fill this in
    end

end
