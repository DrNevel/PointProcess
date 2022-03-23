function [L] = Lambda_eval(xt, wt, thetap, k)
%MAXI_LOGLIKE: maximize the joint log-likelihood of the POINT PROCESS (CIF) to derive the
%optimal parameters for the IG probability distribution.
%   xn: vector of linear system for AR LS estimate
%   wn: matrix of linear system for AR LS estimate
%   thetap: previous uncensored AR parameters estimate
%   k: previous uncensored shape factor estimate

% TODO: definire criterio convergenza
warning('off');

%% [1] Generate handles 
% mu(thetap, xt') RIGHT CENSORING
% mu(thetap, xn) UNCENSORED 
mu = @(thetap,x) x*thetap';

%%%%%%%% Inverse Gaussian pdf
% InvGauss = ...
% @(w,mu,k)          sqrt(k./(2*pi*(w.^3))).* exp(-(k*(w-mu).^2)./(2*w.*(mu.^2))) ;

% Bonvini
InvGauss = ...
@(w,mu,k)          exp(0.5*log(k ./ (2*pi*(w.^3))) - ...
                   0.5*((k*(w - mu).^2) ./ ((mu.^2).*w)));

%%%%%%%% Cumulative distribution of Inverse Gaussian
CumIG= ...
@(w,k,mu)          normcdf(sqrt(k./w).*((w./mu)-1)) + ...
                   exp((2*k./mu) + log(normcdf(-(sqrt(k./w)).*((w./mu)+1))));                       
 
L= InvGauss(wt, mu(thetap, xt'), k) ./ (1-CumIG(wt,k,mu(thetap,xt'))) ;

end