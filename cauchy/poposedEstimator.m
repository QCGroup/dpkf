function [mu_new, P_new] =  poposedEstimator(mu_old, P_old, gamma, C, A, y_t, Q)


% Time update
mu_tmp = A*mu_old;
P_tmp = A*P_old*A' + Q;

% Compute gradient. Approximation is around m_{t | t-1}
v_bar = y_t - C*mu_tmp;
ve = 0;
gr = 2*v_bar / (v_bar^2 + gamma^2);

% New method
Hr = gr/ (v_bar - ve);
P_new = inv(inv(P_tmp) + C'*Hr*C);
mu_new = mu_tmp + P_new*C'*Hr*(y_t - C*mu_tmp - ve);

