function [thetap, k, Loglikel, L] = maxi_loglike_CIF(xn, wn, eta, xt, wt, thetap, k)
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

InvGauss = ...
@(w,mu,k)          exp(0.5*log(k ./ (2*pi*(w.^3))) - ...
                   0.5*((k*(w - mu).^2) ./ ((mu.^2).*w)));

%%%%%%%% Cumulative distribution of Inverse Gaussian
CumIG= ...
@(w,k,mu)          normcdf(sqrt(k./w).*((w./mu)-1)) + ...
                   exp((2*k./mu) + log(normcdf(-(sqrt(k./w)).*((w./mu)+1))));                       
 
%%%%%%%% LogLikelihood
% loglikel = ...to do
% @(w,k,wn,mu)        

%%%%%%%% 1st Order Derivatives of Inverse Gaussian pdf
dkIG = ...
@(k,mu,w)        ( sqrt(1./(8*pi*k*w.^3)) .* exp(-(k*(w-mu).^2)./(2*w.*(mu.^2))) + ...
                  sqrt(k./(2*pi*(w.^3))) .* exp(-(k*(w-mu).^2)./(2*w.*(mu.^2))) .* (-((w-mu).^2)./(2*w.*mu.^2)) );

dthIG = ...         
@(k,mu,w,xn)     ( sqrt(k./(2*pi*(w.^3))) .* exp(-(k*(w-mu).^2)./(2*w.*(mu.^2))) .* ((k*(w-mu))./(w.*mu.^3))).*xn;

%%%%%%%% 1st Order Derivatives of log(Inverse Gaussian)
dkLogIG = ...
@(k,mu,w)        (1/(2*k)) - ( ((w-mu).^2) ./ (2*(mu.^2).*w) );

dthLogIG = ...
@(k,mu,w,xn)     (k * ((w-mu) ./ (mu.^3)) .* xn);


%%%%%%%% 1st Order Derivatives of Cumulative(Inverse Gaussian)
dkCumIG = ...
@(k,mu,w)        ( (normpdf(sqrt(k./w) .* ((w./mu) - 1))) .* ( (1 ./ (2*w.*sqrt(k.*w))) .* ((w./mu) - 1)) +...
                 (2./mu) .* exp((2*k./mu) + log(normcdf(-sqrt(k./w) .* ((w./mu) + 1)))) + ...
                 exp((2*k./mu) + log(normpdf(-sqrt(k./w) .* ((w./mu) + 1)))) .* ( -(1 ./ (2*w.*sqrt(k.*w))) .* ((w./mu) + 1)) ) ;

                        
dthCumIG = ...
@(k,mu,w,xn)     ( (normpdf(sqrt(k./w) .* ((w./mu) - 1))) .* (-sqrt(k./w) .* (w ./ mu.^2)) + ...
                 ((-2*k)./mu.^2) .* exp((2*k./mu) + log(normcdf(-sqrt(k./w) .* ((w./mu) + 1)))) + ...
                 (sqrt(k./w) .* (w./mu.^2)) .*  exp((2*k./mu) + log((normpdf(-sqrt(k./w) .* ((w./mu) + 1)))))   ).*xn;

%% Maximization
        
if ~exist('wt','var')  % UNCENSORED ESTIMATE
        k0 = 1000;
        
        if det(xn'*xn) == 0
            AR0 = pinv(xn'*xn)*(xn'*wn);
        else
            AR0 = inv(xn'*xn)*(xn'*wn);
        end
        
        update = [AR0' k0];
        step = 1;

        while(step<10)     
            % Loglikel before updating
            Loglikel_pre= sum (eta.*( log( InvGauss(wn, mu(update(1:end-1), xn), update(end)) ./ (1-CumIG(wn, update(end), mu(update(1:end-1),xn))) ) )) ;

            % Gradient
            GRAD = eval_grad(CumIG, dkLogIG, dkCumIG, dthLogIG, dthCumIG, mu, update, wn, xn, eta);
            
            % Hessian
            HESS= eval_hess(CumIG, dkLogIG, dkCumIG, dthLogIG, dthCumIG, mu, update, wn, xn, eta);
                       
            FIX = GRAD/HESS;
            update = update - FIX;
            
            % Loglikel after updating
            Loglikel= sum ( log( InvGauss(wn, mu(update(1:end-1), xn), update(end)) ./ (1-CumIG(wn ,update(end), mu(update(1:end-1),xn))) ) ) ; % messo a 0 temporaneamnete per non avere problemi con gli output
            
            if (abs(Loglikel_pre-Loglikel)<10e-5) % exit while loop if loglikel increase is less than 10e-5
               break; 
            end
                        
            step = step + 1;
        end
        
        % Results
        thetap  = update(1:end-1);
        k       = update(end);      
        
else   % RIGHT CENSORING
    
        k0= k;
        update = [thetap , k0];
        step=1;
        WT=[0.005:0.005:wt];
        
        while(step<10)  
            % Loglikel before updating
            Loglikel_pre = sum (eta.*( log( InvGauss(wn, mu(update(1:end-1), xn), update(end)) ./ (1-CumIG(wn ,update(end), mu(update(1:end-1),xn))) ) )) - ...
                    sum ( (InvGauss(WT, mu(update(1:end-1), xt'), update(end)) ./ (1-CumIG(WT,update(end),mu(update(1:end-1),xt'))) ) )*0.005 ;

%            Loglikel_pre = sum ( log( InvGauss(wn, mu(update(1:end-1), xn), update(end)) ./ (1-CumIG(wn ,update(end), mu(update(1:end-1),xn))) ) ) - ...
%                      ( (InvGauss(wt, mu(update(1:end-1), xt'), update(end)) ./ (1-CumIG(wt,update(end),mu(update(1:end-1),xt'))) ) ) * 0.005 ;

            GRAD_u = eval_grad(CumIG, dkLogIG, dkCumIG, dthLogIG, dthCumIG, mu, update, wn, xn,eta);
            GRAD_rc = eval_grad_rc(CumIG, dkIG, dkCumIG, dthIG, dthCumIG, InvGauss, mu, update, wt, xt);
            
            GRAD_rc_tot=[];
            for i=1:length(wn)
                GRAD_temp = eval_grad_rc(CumIG, dkIG, dkCumIG, dthIG, dthCumIG, InvGauss, mu, update, wn(i), xt);
                GRAD_rc_tot = [ GRAD_rc_tot , GRAD_temp] ;
            end
            GRAD_rc_tot = GRAD_rc_tot';
            GRAD_rc_tot = sum(eta.*GRAD_rc_tot);
            
%             GRAD = GRAD_u - GRAD_rc;
            GRAD = GRAD_u - (GRAD_rc_tot + GRAD_rc');
            HESS= eval_hess(CumIG, dkLogIG, dkCumIG, dthLogIG, dthCumIG, mu, update, wn, xn, eta);
           
            FIX = GRAD/HESS;
            
            % temporary update
            updateTemp = update - FIX;
            
            if updateTemp(end)<0 % exit while loop if k parameter is negative
                break;
            end
            
             % Loglikel after updating
            Loglikel= sum ( log( InvGauss(wn, mu(updateTemp(1:end-1), xn), updateTemp(end)) ./ (1-CumIG(wn ,updateTemp(end), mu(updateTemp(1:end-1),xn))) ) ) - ...
                    sum (InvGauss(WT, mu(updateTemp(1:end-1), xt'), updateTemp(end)) ./ (1-CumIG(WT,updateTemp(end),mu(updateTemp(1:end-1),xt'))) )*0.005 ; 
             
            % Loglikel= sum ( log( InvGauss(wn, mu(updateTemp(1:end-1), xn), updateTemp(end)) ./ (1-CumIG(wn ,updateTemp(end), mu(updateTemp(1:end-1),xn))) ) ) - ...
            %                      (InvGauss(wt, mu(updateTemp(1:end-1), xt'), updateTemp(end)) ./ (1-CumIG(wt,updateTemp(end),mu(updateTemp(1:end-1),xt'))) ) * 0.005 ; 
            % 

            if ((Loglikel-Loglikel_pre)>10e-10) 
                update=updateTemp; % update estimate if loglikelihood increase is more than 10e-10
            else
                break; % exit while loop if loglikel increase is less than 10e-10
            end
                
            step = step + 1;
        end
        
        thetap  = update(1:end-1);
        k       = update(end);      
        L= InvGauss(wt, mu(thetap, xt'), k) ./ (sum(1-CumIG(WT,k,mu(thetap,xt')))*0.005) ;
        Loglikel= sum (eta.*( log( InvGauss(wn, mu(thetap, xn), k) ./ (1-CumIG(wn ,k, mu(thetap,xn))) ) )) - ...
                   sum (InvGauss(WT, mu(thetap, xt'), k) ./ (1-CumIG(WT,k,mu(thetap,xt'))) )* 0.005 ; 
%         Loglikel= sum ( log( InvGauss(wn, mu(thetap, xn), k) ./ (1-CumIG(wn ,k, mu(thetap,xn))) ) ) - ...
%                    ((InvGauss(wt, mu(thetap, xt'), k) ./ (1-CumIG(wt,k,mu(thetap,xt'))) ))* 0.005 ; 


end

end