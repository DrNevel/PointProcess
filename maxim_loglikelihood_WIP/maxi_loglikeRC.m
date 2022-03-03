function [thetap, k, Loglikel, L] = maxi_loglikeRC(xn, wn, eta, thetap, k, xt, wt)
%%%%%%%%%% uncensored part funziona, rc no

%MAXI_LOGLIKE: maximize the uncensored joint log-likelihood to derive the
%optimal parameters for the IG probability distribution.
%   xn: vector of linear system for AR LS estimate
%   wn: matrix of linear system for AR LS estimate
%   eta: weighting vector for likelihood moving window
%   thetap: previous uncensored AR parameters estimate
%   k: previous uncensored shape factor estimate
%   xt: ...
%   wt: ...

% TODO: definire criterio convergenza
warning('off');

% if thetap does not exist
%    --> uncensored estimate (init parameters with LS estimate = AR0)
%    --> output thetap, k, steps
% else (thetap exists)
%    --> right censoring estimate (init parameters = thetap)
%    --> output thetap, k, steps, L, loglike
% end

%% [1] Generate handles for the gradient and hessian of the IG joint loglikelihood:
% mu is a vector of estimates conditioned to the history up to u_(i-1)
% mu(thetap, xt') RIGHT CENSORING
% mu(thetap, xn) UNCENSORED 
mu = @(thetap,x) x*thetap';

%                   Inverse Gaussian (noi)
IG = ...
@(wt,mu,k)          sqrt(k./(2*pi*(wt.^3))).* exp(-(k*(wt-mu).^2)./(2*wt.*(mu.^2))) ;

% bonvini
%IG_bonv = ...
%@(wt,mu,k)          exp(0.5*log(k ./ (2*pi*(wt.^3))) - ...
%                    0.5*((k*(wt - mu).^2) ./ ((mu.^2)*wt)));

%                   Cumulative distribution of IG
CIG= ...
@(wt,k,mu)          normcdf(sqrt(k./wt).*((wt./mu)-1)) + ...
                    exp((2*k./mu) + log(normcdf(-(sqrt(k./wt)).*((wt./mu)+1))));    
                
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


%                           Right Censored Likelihood
loglikel_rc = ...
@(w,CIG,loglikel)           loglikel + w(end) * log(1-CIG);

%                           1st Order Derivatives of Right Censored Loglikelihood
dkL_rc = ...
@(w,k,mu,wt,CIG,dkL)        dkL + (... 
                            ((w(end)  / (1 - CIG)) .* ( ...
                            (((wt/mu) -1) / (2*wt*sqrt(k/wt)))* ...
                            (normpdf(sqrt(k/wt) * ((wt/mu) - 1))) + ...
                            ((2/mu) * exp((2*k/mu) + log(normcdf(-sqrt(k/wt) * ((wt/mu) + 1)))))...
                            - (((wt/mu) + 1)/(2*wt*sqrt(k/wt))) * ...
                            (exp((2*k/mu) + log(normpdf(-sqrt(k/wt) * ((wt/mu) + 1)))))...
                            )));
dthL_rc = ...
@(w,k,mu,wt,xt,CIG,dthL)    dthL + (...
                            (w(end)  / (1 - CIG)) .* ( ...
                            (-sqrt(k/wt) * (wt / mu.^2)) * (normpdf(sqrt(k/wt) * ((wt/mu) - 1))) - ...
                            ((2*k)/mu.^2) * (exp((2*k/mu) + log(normcdf(-sqrt(k/wt) * ((wt/mu) + 1))))) + ...
                            (sqrt(k/wt) * (wt/mu.^2)) * exp((2*k/mu) + log(normpdf(-sqrt(k/wt) * ((wt/mu) + 1))))) .* ...
                            xt)';
                         
%% [2] Initialize the parameters for best convergence of Newton-Raphson scheme:
% 
    if ~exist('thetap','var')

        %%%%%%%%%%%%%%%% UNCENSORED ESTIMATE %%%%%%%%%%%%%%%%%
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

    else

        %%%%%%%%%%%%%%%% RIGHT CENSORING ESTIMATE %%%%%%%%%%%%%%%%%
        
        % init
        k0= k;
        update = [thetap k0]; % update(end) is k /// update (1:end-1) is thetap
        step=1;
        h=eps*10000; 

        while(step<10) 
            
            GRAD = [dthL_rc(eta,update(end),mu(update(1:end-1),xt'),wt,xt,CIG(wt,update(end),mu(update(1:end-1),xt')),dthL(eta,update(end),mu(update(1:end-1), xn),wn,xn)) ...
                    dkL_rc(eta,update(end),mu(update(1:end-1),xt'),wt,CIG(wt,update(end),mu(update(1:end-1), xt')),dkL(eta,update(end),mu(update(1:end-1),xn),wn))];
            
            HESS=zeros(length(update));
            for i=1:length(update)
                e=zeros(length(update),1)';
                e(i)=h;
                update_incr=update+e;
                update_decr=update-e;
                              
                GRAD_incr = [dthL_rc(eta,update_incr(end),mu(update_incr(1:end-1),xt'),wt,xt,CIG(wt,update_incr(end),mu(update_incr(1:end-1),xt')),dthL(eta,update_incr(end),mu(update_incr(1:end-1), xn),wn,xn)) ...
                            dkL_rc(eta,update_incr(end),mu(update_incr(1:end-1),xt'),wt,CIG(wt,update_incr(end),mu(update_incr(1:end-1), xt')),dkL(eta,update_incr(end),mu(update_incr(1:end-1),xn),wn))];
                GRAD_decr = [dthL_rc(eta,update_decr(end),mu(update_decr(1:end-1),xt'),wt,xt,CIG(wt,update_decr(end),mu(update_decr(1:end-1),xt')),dthL(eta,update_decr(end),mu(update_decr(1:end-1), xn),wn,xn)) ...
                            dkL_rc(eta,update_decr(end),mu(update_decr(1:end-1),xt'),wt,CIG(wt,update_decr(end),mu(update_decr(1:end-1), xt')),dkL(eta,update_decr(end),mu(update_decr(1:end-1),xn),wn))];        
                hess_i=(GRAD_incr-GRAD_decr)/(2*h);
                HESS(:,i)= hess_i';
            end           
            FIX = GRAD/HESS;
            update = update - FIX;
            step = step + 1;
        end
        % Results
        thetap  = update(1:end-1);
        k       = update(end);
        Loglikel = loglikel_rc(eta, CIG(wt,k,mu(thetap,xt')), loglikel(eta, k, wn, mu(thetap,xn)));
        % Hazard function Lambda
        L = IG(wt,mu(thetap,xt'),k) / (1-CIG(wt,k,mu(thetap,xt')));
    end

end