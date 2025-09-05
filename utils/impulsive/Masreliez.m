function [x_est, P_est] = Masreliez(A, C, Q, y_t, alphaV_1, muV_1, varV_1, muV_2, varV_2, x_prev, P_prev, muW)

%% prediction
mu_tmp = A * x_prev + muW;
P_tmp = A * P_prev * A' + Q;   

%% compute g and G (scalars because measurement is 1-D)
[g, G] = compute_gG(y_t, C, mu_tmp, P_tmp, alphaV_1, muV_1, varV_1, muV_2, varV_2);

%% Masreliez update
x_est = mu_tmp + P_tmp * C' * g;               
P_est = P_tmp - P_tmp * C' * (G * C * P_tmp);           

end
















