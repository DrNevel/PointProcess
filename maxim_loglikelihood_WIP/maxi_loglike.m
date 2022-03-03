function [thetap, k, loglike] = maxi_loglike(xn, wn, eta)
%MAXI_LOGLIKE: maximize the uncensored joint log-likelihood to derive
%   the optimal parameters for the IG probability distribution.
%   xn: vector of linear system for AR LS estimate
%   wn: matrix of linear system for AR LS estimate
%   k_guess: scalar initialization for shape factor
%   eta: weighting vector for likelihood moving window

% TODO: definire criterio convergenza
warning('off');
%% [1] Initialize the parameters for best convergence of Newton-Raphson scheme:
k0 = 1000;

% this kind of check on matrix singularity works well only if AR order
% is not very big:
if det(xn'*xn) == 0
    AR0 = pinv(xn'*xn)*(xn'*wn);
else
    AR0 = inv(xn'*xn)*(xn'*wn);
end

%% [2] Generate handles for the gradient and hessian of the IG joint loglikelihood:

loglike = ...
@(w,k,wn,mu) sum(w.*((1/2)*(log(k)-log(2*pi)-3*log(wn))-(k*(wn-mu).^2)./(2*(mu.^2).*wn)));

% mu is a vector of estimates conditioned to the history up to u_(i-1)
% NOT A SCALAR
mu = @(thetap,xn) xn*thetap' ;

dkL = ...
@(w,k,mu,wn)        sum(w.*((1/(2*k))-(wn-mu).^2./(2*(mu.^2).*wn)));

dthL = ...
@(w,k,mu,wn,xn)     sum(w.*(k*((wn-mu)./mu.^3)'*xn));

d2kL = ...
@(w,k)             -sum(w.*(1/(2*k^2)));

d2thL = ...
@(w,k,mu,wn,xn)    -sum(w.*(k*((3*wn-2*mu)./(mu.^4))))*(xn'*xn);

dkthL = ...
@(w,k,mu,wn,xn)     dthL(w,k,mu,wn,xn)/k;

%% [3] Newton-Raphson set-up:
% GRAD = [dthL(w,k,mu,wn,xn) dkL(w,k,mu,wn)];
% 
% HESS = [d2thL(w,k,mu,wn,xn) dkthL(w,k,mu,wn,xn)'; ...
%         dkthL(w,k,mu,wn,xn) d2kL(w,k)           ];

% Initialize the theta coefficients vector:
update = [AR0' k0];
step = 1;

while(step<10)
    GRAD = [dthL(eta,update(end),mu(update(1:end-1),xn),wn,xn) dkL(eta,update(end),mu(update(1:end-1),xn),wn)];
    HESS = [d2thL(eta,update(end),mu(update(1:end-1),xn),wn,xn) dkthL(eta,update(end),mu(update(1:end-1),xn),wn,xn)'; ...
            dkthL(eta,update(end),mu(update(1:end-1),xn),wn,xn) d2kL(eta,update(end))                              ];
    FIX = GRAD/HESS;
    update = update - FIX;
    step = step + 1;
end

thetap  = update(1:end-1);
k       = update(end);
loglike = loglike(eta, k, wn, mu(thetap,xn));
end