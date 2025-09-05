function [g, G] = compute_gG(y_t, C, mu_tmp, P_tmp, mu, c)
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
    mask = (U > mu);
    U = U(mask);  
    prior = prior(mask);
    
    % Evaluate noise pdf and derivatives (user fills these in)
    v0 = p(U, mu, c);    % p_v(u)
    v1 = dp(U, mu, c);   % derivative wrt u
    v2 = ddp(U, mu, c);  % second derivative wrt u
    
    % Weighted sums (quadrature)
    P  = sum(v0(:) .* prior) * w;
    P1 = sum(v1(:) .* prior) * w;
    P2 = sum(v2(:) .* prior) * w;
    
    % Compute g and G
    g = - P1 / P;
    G = - (P2 / P - (P1 / P)^2);

    % --- Helper functions to be filled by user ---
    function y = p(u, mu, c)
        % Lévy(μ, c) pdf; support u > mu
        y = zeros(size(u));
        m = (u > mu);
        x = u(m) - mu;
        A = sqrt(c/(2*pi));
        y(m) = A .* x.^(-3/2) .* exp(-c./(2*x));
    end
    
    function y = dp(u, mu, c)
        % First derivative w.r.t. u
        y = zeros(size(u));
        m = (u > mu);
        x = u(m) - mu;
        A = sqrt(c/(2*pi));
        f = A .* x.^(-3/2) .* exp(-c./(2*x));
        l1 = -3./(2*x) + c./(2*x.^2);
        y(m) = f .* l1;
    end
    
    function y = ddp(u, mu, c)
        % Second derivative w.r.t. u
        y = zeros(size(u));
        m = (u > mu);
        x = u(m) - mu;
        A = sqrt(c/(2*pi));
        f  = A .* x.^(-3/2) .* exp(-c./(2*x));
        l1 = -3./(2*x) + c./(2*x.^2);
        l2 =  3./(2*x.^2) - c./(x.^3);
        y(m) = f .* (l1.^2 + l2);
    end


end
