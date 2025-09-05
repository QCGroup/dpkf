function r = skew_normal_rnd(xi, w, alpha, N)
% Generate N random numbers from a skew-normal distribution
% with location xi, scale w, and shape alpha

    % Step 1: Generate standard normal variables
    u0 = randn(N, 1);
    v = randn(N, 1);

    % Step 2: Create the skew-normal variates
    delta = alpha / sqrt(1 + alpha^2);
    u1 = delta * u0 + sqrt(1 - delta^2) * v;

    % Step 3: Apply transformation
    z = u0 .* (u1 >= 0) - u0 .* (u1 < 0);

    % Step 4: Scale and shift
    r = xi + w * z;
end
