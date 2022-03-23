function [thetap, k, Loglikel] = maxi_loglike_rc_old(xn, wn, eta)
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
@(w,CIG,loglikel)           loglikel + w * log(1-CIG);

%                           1st Order Derivatives of Right Censored Loglikelihood
%problema: non viene scalare --> 
%(w(end-(length(thetap)-1):end)  / (1 - CIG))      10x1
% altro termine moltiplicativo           1x1
%ho aggiunto sum
dkL_rc = ...
@(w,k,mu,wt,CIG,dkL)        dkL + (... 
                            sum((w(end)  / (1 - CIG)) .* ( ...
                            (((wt/mu) -1) / (2*wt*sqrt(k/wt)))* ...
                            (normpdf(sqrt(k/wt) * ((wt/mu) - 1))) + ...
                            ((2/mu) * exp((2*k/mu) + log(normcdf(-sqrt(k/wt) * ((wt/mu) + 1)))))...
                            - (((wt/mu) + 1)/(2*wt*sqrt(k/wt))) * ...
                            (exp((2*k/mu) + log(normpdf(-sqrt(k/wt) * ((wt/mu) + 1)))))...
                            )));
%%%PROBLEMA DI DIMENSIONI, w è lungo 152x1    
% eta(end-P:end) ???
dthL_rc = ...
@(w,k,mu,wt,xt,CIG,dthL)    dthL + (...
                            (w(end)  / (1 - CIG)) .* ( ...
                            (-sqrt(k/wt) * (wt / mu.^2)) * (normpdf(sqrt(k/wt) * ((wt/mu) - 1))) - ...
                            ((2*k)/mu.^2) * (exp((2*k/mu) + log(normcdf(-sqrt(k/wt) * ((wt/mu) + 1))))) + ...
                            (sqrt(k/wt) * (wt/mu.^2)) * exp((2*k/mu) + log(normpdf(-sqrt(k/wt) * ((wt/mu) + 1))))) .* ...
                            xt)';
                         
%% [2] Initialize the parameters for best convergence of Newton-Raphson scheme:
% 
%     if ~exist('thetap','var')

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

%     else
% 
%         %%%%%%%%%%%%%%%% RIGHT CENSORING ESTIMATE %%%%%%%%%%%%%%%%%
%         
%         % init
%         k0= k;
%         update = [thetap k0];
%         step=1;
%         h=eps;
% 
%         while(step<10) 
%             GRAD = [dthL_rc(eta,update(end),mu(update(1:end-1),xt'),wt,xt,CIG(wt,update(end),mu(update(1:end-1),xt')),dthL(eta,update(end),mu(update(1:end-1), xt'),wt,xt')) ...
%                     dkL_rc(eta,update(end),mu(update(1:end-1),xt'),wt,CIG(wt,update(end),mu(update(1:end-1), xt')),dkL(eta,update(end),mu(update(1:end-1),xt'),wt))];
%             
%             HESS=zeros(length(update));
%             for i=1:length(update)
%                 e=zeros(length(update),1)';
%                 e(i)=h;
%                 update_incr=update+e;
%                 update_decr=update-e;
%                 GRAD_incr = [dthL_rc(eta,update_incr(end),mu(update_incr(1:end-1),xt'),wt,xt,CIG(wt,update_incr(end),mu(update_incr(1:end-1),xt')),dthL(eta,update_incr(end),mu(update_incr(1:end-1), xt'),wt,xt')) ...
%                             dkL_rc(eta,update_incr(end),mu(update_incr(1:end-1),xt'),wt,CIG(wt,update_incr(end),mu(update_incr(1:end-1), xt')),dkL(eta,update_incr(end),mu(update_incr(1:end-1),xt'),wt))];
%                 GRAD_decr = [dthL_rc(eta,update_decr(end),mu(update_decr(1:end-1),xt'),wt,xt,CIG(wt,update_decr(end),mu(update_decr(1:end-1),xt')),dthL(eta,update_decr(end),mu(update_decr(1:end-1), xt'),wt,xt')) ...
%                             dkL_rc(eta,update_decr(end),mu(update_decr(1:end-1),xt'),wt,CIG(wt,update_decr(end),mu(update_decr(1:end-1), xt')),dkL(eta,update_decr(end),mu(update_decr(1:end-1),xt'),wt))];        
%                 hess_i=(GRAD_incr-GRAD_decr)/(2*h);
%                 HESS(:,i)= hess_i';
%             end
% %             dkL_rc_xincr= dkL_rc(eta ,k ,mu(thetap+h, xt') ,wt ,CIG(wt, update(end), mu(thetap+h, xt')) ,dkL(eta, k, mu(thetap+h, xt'), wn));
% %             dkL_rc_x=dkL_rc(eta ,k ,mu(thetap, xt') ,wt ,CIG(wt, update(end), mu(thetap, xt')) ,dkL(eta, k, mu(thetap, xt'), wn));
% %           
% %             dkL_rc_xincr= dkL_rc(eta ,update(end)+h ,mu(update(1:end-1)+h, xt') ,wt ,CIG(wt, update(end)+h, mu(update(1:end-1), xt')) ,dkL(eta, update(end)+k, mu(update(1:end-1), xt'), wn));
% %             dkL_rc_x=dkL_rc(eta ,update(end) ,mu(update(1:end-1), xt') ,wt ,CIG(wt, update(end), mu(update(1:end-1), xt')) ,dkL(eta, update(end), mu(update(1:end-1), xt'), wn));
% %             
% %             dthL_rc_xincr=dthL_rc(eta,update(end),mu(update(1:end-1),xt'),wt,xt,CIG(wt,update(end),mu(update(1:end-1),xt')),dthL(eta,update(end),mu(update(1:end-1), xt'),wt,xt'));
% %             dthL_rc_x=dthL_rc(eta,update(end),mu(update(1:end-1),xt'),wt,xt,CIG(wt,update(end),mu(update(1:end-1),xt')),dthL(eta,update(end),mu(update(1:end-1), xt'),wt,xt'));
% %             
% %             
% %             cross = ((dkL_rc_xincr - dkL_rc_x) / 2*step) + ((dthL_rc_xincr - dthL_rc_x) /2*step);
% %             th2 = ((dthL_rc_xincr - dthL_rc_x) /2*step) + ((dthL_rc_xincr - dthL_rc_x) /2*step);
% %             k2 = ((dkL_rc_xincr - dkL_rc_x) / 2*step) + ((dkL_rc_xincr - dkL_rc_x) / 2*step);
% %             
% %             th2=[];
% %             for i=1:length(xt)
% %                 temp=update(i);
% %                 update(i)=update(i)+h;
% %                 th2=[th2 ; dthL_rc(eta,update(end),mu(update(1:end-1),xt'),wt,xt,CIG(wt,update(end),mu(update(1:end-1),xt')),dthL(eta,update(end),mu(update(1:end-1), xt'),wt,xt'))];
% %                 update(i)=temp;
% %             end
% %             
% %             th2=th2 - dthL_rc(eta,update(end),mu(update(1:end-1),xt'),wt,xt,CIG(wt,update(end),mu(update(1:end-1),xt')),dthL(eta,update(end),mu(update(1:end-1), xt'),wt,xt'));
% %             
% %             HESS = [th2 ...
% %                    cross'; ...
% %                    cross...
% %                    k2];
%             
%             %FIX = GRAD/HESS;
%             FIX = GRAD*pinv(HESS);
%             update = update - FIX;
%             step = step + 1;
%         end
% 
%         % Results
%         thetap  = update(1:end-1);
%         k       = update(end);
%         Loglikel = loglikel_rc(eta, CIG(wt,k,mu(thetap,xt')), loglikel(eta, k, wt, mu(thetap,xt')));
%         % Hazard function Lambda
%         L = IG(wt,mu(thetap,xt'),k) / (1-CIG(wt,k,mu(thetap,xt')));
%         
%     end

end