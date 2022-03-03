function [thetap, k, Loglikel] = maxi_loglikeRC_FINAL(xn, wn)
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

%% [1] Generate handles 
% mu(thetap, xt') RIGHT CENSORING
% mu(thetap, xn) UNCENSORED 
mu = @(thetap,x) x*thetap';

%                   Inverse Gaussian probability density function
InvGauss = ...
@(w,mu,k)          sqrt(k./(2*pi*(w.^3))).* exp(-(k*(w-mu).^2)./(2*w.*(mu.^2))) ;

% Bonvini
%IG_bonv = ...
%@(wt,mu,k)          exp(0.5*log(k ./ (2*pi*(wt.^3))) - ...
%                    0.5*((k*(wt - mu).^2) ./ ((mu.^2)*wt)));

%                   Cumulative distribution of Inverse Gaussian
CumIG= ...
@(w,k,mu)          normcdf(sqrt(k./w).*((w./mu)-1)) + ...
                   exp((2*k./mu) + log(normcdf(-(sqrt(k./w)).*((w./mu)+1))));    
                
%                   LogLikelihood
% loglikel = ...to do
% @(w,k,wn,mu)        

%                   1st Order Derivatives of Inverse Gaussian probability
%                   density function
dkIG = ...
@(k,mu,w)        ((1./(4.*(w.^3).*pi.*sqrt(k./(2*pi*w.^3)))) .* exp(-(k*(w-mu).^2)./(2*w.*(mu.^2))) + ...
                  sqrt(k./(2*pi*(w.^3))) .* exp(-(k*(w-mu).^2)./(2*w.*(mu.^2))) .* ((-k*(w-mu).^2)./(2*w.*mu.^2)));

dthIG = ...         
@(k,mu,w,xn)     (sqrt(k./(2*pi*(w.^3))) .* exp(-(k*(w-mu).^2)./(2*w.*(mu.^2))) .* ((k*(w-mu))./(w.*mu.^3)))'*xn;

%                   1st Order Derivatives of log(Inverse Gaussian)
dkLogIG = ...
@(k,mu,w)        ((1/(2*k))-(w-mu).^2./(2*(mu.^2).*w));

dthLogIG = ...
@(k,mu,w,xn)     (k*((w-mu)./mu.^3)'*xn);


%                   1st Order Derivatives of Cumulative(Inverse Gaussian)
dkCumIG = ...
@(k,mu,w)        (normpdf(sqrt(k./w) .* ((w./mu) - 1))) .* ((1 ./ (2*sqrt(k./w).*w)) .* ((w./mu) - 1)) +...
                 (2./mu) .* exp((2*k./mu) + log(normcdf(-sqrt(k./w) .* ((w./mu) + 1)))) + ...
                 exp((2*k./mu) + log(normpdf(-sqrt(k./w) .* ((w./mu) + 1)))) .* ((1 ./ (2*sqrt(k./w).*w)) .* ((w./mu) - 1));

                        
dthCumIG = ...
@(k,mu,w,xn)     ( (normpdf(sqrt(k/w) * ((w/mu) - 1))) .* (-sqrt(k/w) * (w / mu.^2)) + ...
                 ((-2*k)/mu.^2) .* exp((2*k/mu) + log(normcdf(-sqrt(k/w) * ((w/mu) + 1)))) + ...
                 (sqrt(k/w) * (w/mu.^2)) .*  exp((2*k/mu) + log((normpdf(-sqrt(k/w) * ((w/mu) + 1)))))   )*xn;

%% Maximization
        k0 = 1000;
        
        if det(xn'*xn) == 0
            AR0 = pinv(xn'*xn)*(xn'*wn);
        else
            AR0 = inv(xn'*xn)*(xn'*wn);
        end
        
        update = [AR0' k0];
        step = 1;
%         incr=1e-10;

        while(step<10) 
            % integral_dN( d_LogIG + (d_CumIG / (1-CumIG)) )   -   ...
            % integral_dt( (d_IG / (1-CumIG)) + ((IG * CumIG) / ((1-CumIG)^2) )
            
            % Gradiente
            kLogIG= dkLogIG(update(end),mu(update(1:end-1),xn),wn);
            thLogIG= dthLogIG(update(end),mu(update(1:end-1),xn),wn,xn);
            
            kIG=dkIG(update(end),mu(update(1:end-1),xn),wn) ;
            thIG=dthIG(update(end),mu(update(1:end-1),xn),wn, xn);
            
            kCumIG=dkCumIG(update(end),mu(update(1:end-1),xn), wn);
            thCumIG=dthCumIG(update(end),mu(update(1:end-1),xn), wn, xn);
            
            CIG=CumIG(wn,update(end),mu(update(1:end-1),xn));
            IG=InvGauss(wn,mu(update(1:end-1),xn),update(end));
            
            grad_k= sum( kLogIG + ( kCumIG ./ (1-CIG) ) ) - ...
                    sum( (kIG ./ (1-CIG))  + ( (IG .* kCumIG) ./ ((1-CIG).^2) ) ) ;
            grad_th= sum( thLogIG + ( thCumIG ./ (1-CIG) ) ) - ...
                    sum( (thIG ./ (1-CIG))  + ( (IG * thCumIG) ./ ((1-CIG).^2) ) ) ;
                
            GRAD = [ grad_th , grad_k ];
            
            % Hessiana approssimata
            HESS=zeros(length(update));
            for i=1:length(update)
                h_i=zeros(length(update),1)';
                incr_i=nthroot(eps,3)*(1+abs(update(i)));
                h_i(i)=incr_i;
                
                update_i=update+h_i;
                
                kLogIG= dkLogIG(update_i(end),mu(update_i(1:end-1),xn),wn);
                thLogIG=dthLogIG(update_i(end),mu(update_i(1:end-1),xn),wn,xn);
            
                kIG=dkIG(update_i(end),mu(update_i(1:end-1),xn), wn) ;
                thIG=dthIG(update_i(end),mu(update_i(1:end-1),xn), wn, xn);
                
                kCumIG=dkCumIG(update_i(end),mu(update_i(1:end-1),xn), wn);
                thCumIG=dthCumIG(update_i(end),mu(update_i(1:end-1),xn), wn, xn);
                
                CIG=CumIG(wn,update_i(end),mu(update_i(1:end-1),xn));
                IG=InvGauss(wn,mu(update_i(1:end-1),xn),update_i(end));
                                
                grad_k_i= sum( kLogIG + ( kCumIG ./ (1-CIG) ) ) - ...
                             sum( (kIG ./ (1-CIG))  + ( (IG .* kCumIG) ./ ((1-CIG).^2) ) ) ;
                grad_th_i= sum( thLogIG + ( thCumIG ./ (1-CIG) ) ) - ...
                              sum( (thIG ./ (1-CIG))  + ( (IG * thCumIG) ./ ((1-CIG).^2) ) ) ;
                GRAD_i= [ grad_th_i , grad_k_i ];
                
                for j=1:length(update)
                    h_j=zeros(length(update),1)';
                    incr_j=nthroot(eps,3)*(1+abs(update(j)));
                    h_j(j)=incr_j;
                    
                    update_j=update+h_j;
                
                    kLogIG= dkLogIG(update_j(end),mu(update_j(1:end-1),xn),wn);
                    thLogIG=dthLogIG(update_j(end),mu(update_j(1:end-1),xn),wn,xn);

                    kIG=dkIG(update_j(end),mu(update_j(1:end-1),xn), wn) ;
                    thIG=dthIG(update_j(end),mu(update_j(1:end-1),xn), wn, xn);

                    kCumIG=dkCumIG(update_j(end),mu(update_j(1:end-1),xn), wn);
                    thCumIG=dthCumIG(update_j(end),mu(update_j(1:end-1),xn), wn, xn);

                    CIG=CumIG(wn,update_j(end),mu(update_j(1:end-1),xn));
                    IG=InvGauss(wn,mu(update_j(1:end-1),xn),update_j(end));

                    grad_k_j= sum( kLogIG + ( kCumIG ./ (1-CIG) ) ) - ...
                                 sum( (kIG ./ (1-CIG))  + ( (IG .* kCumIG) ./ ((1-CIG).^2) ) ) ;
                    grad_th_j= sum( thLogIG + ( thCumIG ./ (1-CIG) ) ) - ...
                                  sum( (thIG ./ (1-CIG))  + ( (IG * thCumIG) ./ ((1-CIG).^2) ) ) ;
                    GRAD_j= [ grad_th_j , grad_k_j ];
                                    
                    HESS(i,j)=( (GRAD_i(i)-GRAD(i))/(2*incr_i) )+( (GRAD_j(j)-GRAD(j))/(2*incr_j) );
                    
                end
            end
                       
            FIX = GRAD/HESS;
            update = update - FIX;
            step = step + 1;
        end
        
        % Results
        thetap  = update(1:end-1);
        k       = update(end);
        
        % da aggiungere:
        %Loglikel = loglikel(eta, k, wn, mu(thetap,xn));
        %L = IG(wt,mu(thetap,xt'),k) / (1-CIG(wt,k,mu(thetap,xt')));
        Loglikel=0; % messo a 0 temporaneamnete per non avere problemi con gli output


%% Appunti: uncensored (wn, xn)
% wrt theta
% sum(eta*( dthLIG + (dthCIG / (1-CIG)) ));
% dthLIG(k,mu(thetap',xn),wn,xn)
% 1-CIG(wn,k,mu(thetap',xn))
% dthCIG(k,mu(thetap',xn), wn, xn)
% dthIG(k,mu(thetap',xn),wn,xn) 

% U_th=sum((dthLIG(k,mu(thetap',xn),wn,xn) + ( dthCIG(k,mu(thetap',xn), wn, xn) ./ (1-CIG(wn,k,mu(thetap',xn)))) ) );

% wrt k
% dkLIG(k,mu(thetap',xn),wn)
% 1-CIG(wn,k,mu(thetap',xn))
% dkCIG(k,mu(thetap',xn), wn)
% dkIG(k,mu(thetap',xn),wn)
% 
% U_k=sum(( dkLIG(k,mu(thetap',xn),wn)' + ( dkCIG(k,mu(thetap',xn), wn) ./ (1-CIG(wn,k,mu(thetap',xn)))' )    ) );
% 
% GRAD= [ U_th , U_k ];
%% Appunti: censored (wt,xt)
% wrt theta
% sum(eta*( (dthIG*(1-CIG) + IG*(dthCIG)) / (1-CIG)^2 ));
% dthIG(k,mu(thetap',xt'),wt,xt) 
% 1-CIG(wt,k,mu(thetap',xt'))
% IG(wt,mu(thetap',xt'),k)
% dthCIG(k,mu(thetap',xt'), wt, xt)
% 
% RC_th=sum(eta.*( (dthIG(k,mu(thetap',xt'),wt,xt)*(1-CIG(wt,k,mu(thetap',xt'))) + IG(wt,mu(thetap',xt'),k)*(dthCIG(k,mu(thetap',xt'), wt, xt))) / (1-CIG(wt,k,mu(thetap',xt')))^2 ));
% 
% wrt k
%RC_k...

end