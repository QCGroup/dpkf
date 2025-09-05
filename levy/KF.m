function [mu_new, P_new] = KF(mu_old, P_old, muV, y_t, Q, R, A, C)

% Time update
mu_tmp = A*mu_old;
P_tmp = A*P_old*A' + Q;

% Measurement update
P_new = inv( inv(P_tmp) + C'*inv(R)*C );
mu_new = mu_tmp + P_new*C'*inv(R)*(y_t - C*mu_tmp - muV);