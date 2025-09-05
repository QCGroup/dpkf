function [g, G] = compute_gG(y_t, C, mu_tmp, P_tmp, lambda)
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
    mask = (U > 0);
    U = U(mask);           % keep only residuals inside support
    prior = prior(mask);   % keep corresponding Gaussian prior weights
    
    % Evaluate noise pdf and derivatives (user fills these in)
    v0 = p(U, lambda);    % p_v(u)
    v1 = dp(U, lambda);   % derivative wrt u
    v2 = ddp(U, lambda);  % second derivative wrt u
    
    % Weighted sums (quadrature)
    P  = sum(v0(:) .* prior) * w;
    P1 = sum(v1(:) .* prior) * w;
    P2 = sum(v2(:) .* prior) * w;
    
    % Compute g and G
    g = - P1 / P;
    G = - (P2 / P - (P1 / P)^2);

    % --- Helper functions to be filled by user ---
    function y = p(u, lambda)
    % Exponential PDF (rate lambda), support u >= 0
    % p(u) = lambda * exp(-lambda*u) for u >= 0, else 0
    y = lambda*ones(size(u));
    m = (u >= 0);
    y(m) = lambda .* exp(-lambda .* u(m));
    end
    
    function y = dp(u, lambda)
        % First derivative wrt u
        % For u > 0: p'(u) = -lambda^2 * exp(-lambda*u)
        % For u < 0: p'(u) = 0
        % At u = 0: use right-hand derivative = -lambda^2
        y = -lambda.^2*ones(size(u));
        m = (u >= 0);  % include u==0 with right derivative
        y(m) = - (lambda.^2) .* exp(-lambda .* u(m));
    end
    
    function y = ddp(u, lambda)
        % Second derivative wrt u (right-hand)
        % For u > 0: p''(u) = lambda^3 * exp(-lambda*u)
        % For u < 0: p''(u) = 0
        % At u = 0: use right-hand value = lambda^3
        y = lambda.^3*ones(size(u));
        m = (u >= 0);
        y(m) = (lambda.^3) .* exp(-lambda .* u(m));
    end


end
