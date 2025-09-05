function [g, G] = compute_gG(y_t, C, mu_tmp, P_tmp, alphaV_1, muV_1, varV_1, muV_2, varV_2)
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
    v0 = p(U, alphaV_1, muV_1, varV_1, muV_2, varV_2);    % p_v(u)
    v1 = dp(U, alphaV_1, muV_1, varV_1, muV_2, varV_2);   % derivative wrt u
    v2 = ddp(U, alphaV_1, muV_1, varV_1, muV_2, varV_2);  % second derivative wrt u
    
    % Weighted sums (quadrature)
    P  = sum(v0(:) .* prior) * w;
    P1 = sum(v1(:) .* prior) * w;
    P2 = sum(v2(:) .* prior) * w;
    
    % Compute g and G
    g = - P1 / P;
    G = - (P2 / P - (P1 / P)^2);

    % --- Helper functions to be filled by user ---
    function y = p(u, alphaV_1, muV_1, varV_1, muV_2, varV_2)
    % PDF of alpha1*N(mu1,var1) + alpha2*N(mu2,var2)
    a1 = alphaV_1; a2 = 1 - a1;
    s1 = sqrt(varV_1); s2 = sqrt(varV_2);
    phi1 = (1./(sqrt(2*pi)*s1)) .* exp(-0.5*((u - muV_1).^2)./varV_1);
    phi2 = (1./(sqrt(2*pi)*s2)) .* exp(-0.5*((u - muV_2).^2)./varV_2);
    y = a1.*phi1 + a2.*phi2;
    end
    
    function y = dp(u, alphaV_1, muV_1, varV_1, muV_2, varV_2)
        % First derivative wrt u
        a1 = alphaV_1; a2 = 1 - a1;
        s1 = sqrt(varV_1); s2 = sqrt(varV_2);
        phi1 = (1./(sqrt(2*pi)*s1)) .* exp(-0.5*((u - muV_1).^2)./varV_1);
        phi2 = (1./(sqrt(2*pi)*s2)) .* exp(-0.5*((u - muV_2).^2)./varV_2);
        dphi1 = phi1 .* ( -(u - muV_1)./varV_1 );
        dphi2 = phi2 .* ( -(u - muV_2)./varV_2 );
        y = a1.*dphi1 + a2.*dphi2;
    end
    
    function y = ddp(u, alphaV_1, muV_1, varV_1, muV_2, varV_2)
        % Second derivative wrt u
        a1 = alphaV_1; a2 = 1 - a1;
        s1 = sqrt(varV_1); s2 = sqrt(varV_2);
        phi1 = (1./(sqrt(2*pi)*s1)) .* exp(-0.5*((u - muV_1).^2)./varV_1);
        phi2 = (1./(sqrt(2*pi)*s2)) .* exp(-0.5*((u - muV_2).^2)./varV_2);
        d2phi1 = phi1 .* ( ((u - muV_1).^2)./(varV_1.^2) - 1./varV_1 );
        d2phi2 = phi2 .* ( ((u - muV_2).^2)./(varV_2.^2) - 1./varV_2 );
        y = a1.*d2phi1 + a2.*d2phi2;
    end
end
