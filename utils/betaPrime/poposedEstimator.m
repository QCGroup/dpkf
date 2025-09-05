function [mu_new, P_new] =  poposedEstimator(mu_old, P_old, alpha, beta, C, A, y_t, Q)

% Time update
mu_tmp = A*mu_old;
P_tmp = A*P_old*A' + Q;

% Compute gradient. Approximation is around m_{t | t-1}
v_bar = y_t - C*mu_tmp;
if v_bar < 0 
    v_bar = 1e-8;
end
ve = (alpha-1)/(beta + 1);
gr = -(alpha - 1)/v_bar + (alpha + beta)/(1 + v_bar);

% New method
Hr = gr/ (v_bar - ve);
P_new = inv(inv(P_tmp) + C'*Hr*C);
mu_new = mu_tmp + P_new*C'*gr;

