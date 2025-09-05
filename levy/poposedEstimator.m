function [mu_new, P_new] =  poposedEstimator(mu_old, P_old, mu, c, C, A, y_t, Q)

% Time update
mu_tmp = A*mu_old;
P_tmp = A*P_old*A' + Q;

% Compute gradient. Approximation is around m_{t | t-1}
v_bar = y_t - C*mu_tmp;
if v_bar < mu
    v_bar = mu + 1e-6;
end
ve = mu + c/3;
gr = 3/2/(v_bar - mu) - c/2/(v_bar - mu)^2;

% New method
Hr = gr/ (v_bar - ve);
P_new = inv(inv(P_tmp) + C'*Hr*C);
mu_new = mu_tmp + P_new*C'*Hr*(y_t - C*mu_tmp - ve);

