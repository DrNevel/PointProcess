function [L, Loglikel_rc, InvGauss, CumIG] = estim(eta, thetap, k, xt, wt, wn, xn)
%% [1] Generate handles 
mu = @(thetap,x) x*thetap';

%                   Inverse Gaussian (noi)
IG = ...
@(wt,mu,k)          sqrt(k./(2*pi*wt.^3)).* exp(-(k*(wt-mu).^2)./(2*wt.*(mu.^2))) ;

% bonvini
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

%                           Right Censored Likelihood
loglikel_rc = ...
@(w,CIG,loglikel)           loglikel + w(end) * log(1-CIG);
                         
%% Calcola loglikel right censoring, inverse gaussian, cif of ig, lambda (hazard funct)

        Loglikel_rc = loglikel_rc(eta, CIG(wt,k,mu(thetap,xt')), loglikel(eta, k, wn, mu(thetap,xn)));
        InvGauss = IG(wt,mu(thetap,xt'),k);
        CumIG = CIG(wt,k,mu(thetap,xt'));
        L = IG(wt,mu(thetap,xt'),k) / (1-CIG(wt,k,mu(thetap,xt')));
    
end