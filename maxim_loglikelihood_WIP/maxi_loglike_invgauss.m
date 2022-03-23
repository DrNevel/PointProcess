function [thetap, k, Loglikel] = maxi_loglike_invgauss(xn, wn, eta)
%% [1] Generate handles for the gradient and hessian of the IG joint loglikelihood:
% mu is a vector of estimates conditioned to the history up to u_(i-1)
% mu(thetap, xt') RIGHT CENSORING
% mu(thetap, xn) UNCENSORED 
mu = @(thetap,x) x*thetap';

%                   Inverse Gaussian (noi)
IG = ...
@(wt,mu,k)          sqrt(k./(2*pi*wt.^3)).* exp(-(k*(wt-mu).^2)./(2*wt.*(mu.^2))) ;

%IG_bonv = ...
%@(wt,mu,k)          exp(0.5*log(k ./ (2*pi*wt.^3)) - ...
%                    0.5*(k*(wt - mu).^2) ./ ((mu.^2)*wt));

%                   Cumulative distribution of IG
CIG= ...
@(wt,k,mu)          normcdf(sqrt(k./wt).*((wt./mu)-1)) + ...
                    exp(2*k./mu + log(normcdf(-(sqrt(k./wt)).*((wt./mu)+1))));    
                
%                   Uncensored Likelihood
loglikel = ...
@(w,k,wn,mu)        sum(w.*((1/2)*(log(k)-log(2*pi)-3*log(wn))-(k*(wn-mu).^2)./(2*(mu.^2).*wn)));

%                   1st and 2nd Order Derivatives of Uncensored Loglikelihood
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

% init
        k0 = 1000;
        % this kind of check on matrix singularity works well only if AR order
        % is not very big:
        if det(xn'*xn) == 0
            AR0 = pinv(xn'*xn)*(xn'*wn);
        else
            AR0 = inv(xn'*xn)*(xn'*wn);
        end

        %%%% NR optimization
        % GRAD = [dthL(w,k,mu,wn,xn) dkL(w,k,mu,wn)];
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

        % Results
        thetap  = update(1:end-1);
        k       = update(end);
        Loglikel = loglikel(eta, k, wn, mu(thetap,xn));
end