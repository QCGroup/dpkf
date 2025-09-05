function [mu_new, P_new] =  poposedEstimator(mu_old, P_old, alphaV_1, muV_1, varV_1, alphaV_2, muV_2, varV_2, C, A, y_t, Q, ve)


% Time update
mu_tmp = A*mu_old;
P_tmp = A*P_old*A' + Q;

% Compute gradient. Approximation is around m_{t | t-1}
v_bar = y_t - C*mu_tmp;
gr = gradr(alphaV_1,alphaV_2,muV_1,muV_2,v_bar,varV_1,varV_2);

% New method
tmp1 = gr/ (v_bar - ve(1));
tmp2 = gr/ (v_bar - ve(2));
if tmp1*tmp2 > 0 
    Hr = max(tmp1, tmp2);
else
    tmp = [tmp1; tmp2];
    Hr = tmp(tmp > 0);
end
P_new = inv(inv(P_tmp) + C'*Hr*C);
mu_new = mu_tmp + P_new*C'*gr;

%% Gradient

function out = gradr(alphaV_1,alphaV_2,muV_1,muV_2,vbar,varV_1,varV_2)
logphi1 = -0.5*log(2*pi*varV_1) - (vbar - muV_1)^2 / (2*varV_1);
    logphi2 = -0.5*log(2*pi*varV_2) - (vbar - muV_2)^2 / (2*varV_2);
    
    logwphi1 = log(alphaV_1) + logphi1;
    logwphi2 = log(alphaV_2) + logphi2;
    
    % log-sum-exp
    maxlog = max(logwphi1, logwphi2);
    logS = maxlog + log( exp(logwphi1 - maxlog) + exp(logwphi2 - maxlog) );
    
    % responsibilities in log-domain
    r1 = exp(logwphi1 - logS);
    r2 = exp(logwphi2 - logS);
    
    % derivative
    out = r1*(vbar - muV_1)/varV_1 + r2*(vbar - muV_2)/varV_2;
end

end