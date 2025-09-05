function [mu_new, P_new] =  poposedEstimator(mu_old, P_old, alpha, scale, C, A, y_t, Q, ve)


% Time update
mu_tmp = A*mu_old;
P_tmp = A*P_old*A' + Q;

% Compute the gradient.
v_bar = y_t - C*mu_tmp;
gr = v_bar/scale^2 - sqrt(2/pi)*alpha/(scale*erfcx(-v_bar*alpha/sqrt(2*scale^2)));

% New method
Hr = gr/ (v_bar - ve);
P_new = inv(inv(P_tmp) + C'*Hr*C);
mu_new = mu_tmp + P_new*C'*gr;