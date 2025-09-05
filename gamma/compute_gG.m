function [g, G] = compute_gG(y_t, C, mu_tmp, P_tmp, alpha, theta)
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
    mask = U > 0;
    U = U(mask);           % keep only residuals inside support
    prior = prior(mask);   % keep corresponding Gaussian prior weights
    
    % Evaluate noise pdf and derivatives (user fills these in)
    v0 = p(U, alpha, theta);    % p_v(u)
    v1 = dp(U, alpha, theta);   % derivative wrt u
    v2 = ddp(U, alpha, theta);  % second derivative wrt u
    
    % Weighted sums (quadrature)
    P  = sum(v0(:) .* prior) * w;
    P1 = sum(v1(:) .* prior) * w;
    P2 = sum(v2(:) .* prior) * w;
    
    % Compute g and G
    g = - P1 / P;
    G = - (P2 / P - (P1 / P)^2);

    % --- Helper functions to be filled by user ---
    function y = p(u, alpha, theta)
    % Gamma(α, θ) PDF: u^(α-1) * exp(-u/θ) / (Γ(α) θ^α),  support u >= 0
    y = zeros(size(u));
    m = (u > 0);
    % stable log form
    logp = (alpha-1).*log(u(m)) - u(m)./theta - alpha.*log(theta) - gammaln(alpha);
    y(m) = exp(logp);
    end
    
    function y = dp(u, alpha, theta)
        % First derivative wrt u:
        % p'(u) = p(u) * ( (α-1)/u - 1/θ ), for u>0; 0 for u<=0
        y = 1e-8*ones(size(u));
        m = (u > 0);
        pu = p(u(m), alpha, theta);
        y(m) = pu .* ( (alpha-1)./u(m) - 1./theta );
    end
    
    function y = ddp(u, alpha, theta)
        % Second derivative wrt u:
        % p''(u) = p(u)*[ ((α-1)/u - 1/θ)^2 - (α-1)/u^2 ]
        %        = p(u)*( (α-1)(α-2)/u^2 - 2(α-1)/(uθ) + 1/θ^2 ), for u>0; 0 else
        y = 1e-8*ones(size(u));
        m = (u > 0);
        pu = p(u(m), alpha, theta);
        y(m) = pu .* ( ((alpha-1).*(alpha-2))./(u(m).^2) ...
                     - 2*(alpha-1)./(u(m).*theta) ...
                     + 1./(theta.^2) );
    end

end
