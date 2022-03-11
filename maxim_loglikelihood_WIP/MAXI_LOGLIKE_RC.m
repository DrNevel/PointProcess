

function [thetap, k, Loglikel] = MAXI_LOGLIKE_RC(xn, wn)
%MAXI_LOGLIKE: maximize the joint log-likelihood of the POINT PROCESS (CIF) to derive the
%optimal parameters for the IG probability distribution.
%   xn: vector of linear system for AR LS estimate
%   wn: matrix of linear system for AR LS estimate
%   thetap: previous uncensored AR parameters estimate
%   k: previous uncensored shape factor estimate

% TODO: definire criterio convergenza
% warning('off');

%% [1] Generate handles 
% mu(thetap, xt') RIGHT CENSORING
% mu(thetap, xn) UNCENSORED 
mu = @(thetap,x) x*thetap';

%%%%%%%% Inverse Gaussian pdf
InvGauss = ...
@(w,mu,k)          sqrt(k./(2*pi*(w.^3))).* exp(-(k*(w-mu).^2)./(2*w.*(mu.^2))) ;

% Bonvini
%IG_bonv = ...
%@(w,mu,k)          exp(0.5*log(k ./ (2*pi*(w.^3))) - ...
%                    0.5*((k*(wt - mu).^2) ./ ((mu.^2)*w)));

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
        k0 = 1000;
        
        if det(xn'*xn) == 0
            AR0 = pinv(xn'*xn)*(xn'*wn);
        else
            AR0 = inv(xn'*xn)*(xn'*wn);
        end
        
        update = [AR0' k0];
        step = 1;

        while(step<10)
                       
            % Gradiente
            GRAD = eval_grad(CumIG, dkLogIG, dkCumIG, dthLogIG, dthCumIG, mu, update, wn, xn);
            
            % Hessiana approssimata
            HESS= eval_hess(CumIG, dkLogIG, dkCumIG, dthLogIG, dthCumIG, mu, update, wn, xn);
                       
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

end

function GRAD = eval_grad(CumIG, dkLogIG, dkCumIG, dthLogIG, dthCumIG, mu, update, wn, xn)
    kLogIG= dkLogIG(update(end),mu(update(1:end-1),xn),wn);
    thLogIG= dthLogIG(update(end),mu(update(1:end-1),xn),wn,xn);
    kCumIG=dkCumIG(update(end),mu(update(1:end-1),xn), wn);
    thCumIG=dthCumIG(update(end),mu(update(1:end-1),xn), wn, xn);
    CIG=CumIG(wn,update(end),mu(update(1:end-1),xn));

%     grad_k= sum( kLogIG + ( kCumIG ./ (1-CIG) ) ) ;
%     grad_th= sum( thLogIG + ( thCumIG ./ (1-CIG) ) );

    %GRAD = [ grad_th , grad_k ];
%     GRAD = [sum(thLogIG) , sum(kLogIG)] + [ sum(thCumIG) , sum(kCumIG) ] ./ sum(1-CIG);
%     GRAD = sum([ thLogIG , kLogIG ] + ( [ thCumIG , kCumIG ] ./ (1-CIG) ));
%     GRAD = sum([ thLogIG , kLogIG ]);
    GRAD = sum([ thLogIG , kLogIG ]) +  sum([ thCumIG , kCumIG ]) ./ sum(1-CIG);
end

function HESS = eval_hess(CumIG, dkLogIG, dkCumIG, dthLogIG, dthCumIG, mu, update, wn, xn)
    
HESS=zeros(length(update)); 
incr=update;
decr=update;

    for i=1:length(update)
%         step=nthroot(eps,3)*(abs(update(i)));
        incr=update;
        decr=update;
        step=abs(update(i))*(nthroot(eps,3));
        
        incr(i) = incr(i) + step;
        decr(i) = decr(i) - step;
        
        GRAD_incr = eval_grad(CumIG, dkLogIG, dkCumIG, dthLogIG, dthCumIG, mu, incr, wn, xn);
        GRAD_decr = eval_grad(CumIG, dkLogIG, dkCumIG, dthLogIG, dthCumIG, mu, decr, wn, xn);
        HESS(i,:)= (GRAD_incr - GRAD_decr ) ./ (2*step);

%         for j=1:length(update)          
%             HESS(i,j)= (GRAD_incr(j) - GRAD_decr(j) ) ./ (2*step) ;
%         end

    end
end