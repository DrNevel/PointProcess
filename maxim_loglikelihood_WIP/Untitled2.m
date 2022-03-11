%% uncensored?
grad_LogIG = [ dthLogIG(update(end),mu(update(1:end-1),xn),wn,xn) , dkLogIG(update(end),mu(update(1:end-1),xn),wn) ];

grad_CUM = [ dthCumIG(update(end),mu(update(1:end-1),xn), wn, xn) dkCumIG(update(end),mu(update(1:end-1),xn), wn) ];

CUM = CumIG(wn,update(end),mu(update(1:end-1),xn));

GRAD_U = (grad_LogIG + grad_CUM ./ (1-CUM)) ;
GRAD_U=sum(GRAD_U .* 0.005)

%% censored?

% grad_IG_rc = [ dkIG(update(end),mu(update(1:end-1),xt'),wt) ,dthIG(update(end),mu(update(1:end-1),xt'),wt, xt') ];
% 
% grad_CUM_rc = [ dkCumIG(update(end),mu(update(1:end-1),xt'), wt) , dthCumIG(update(end),mu(update(1:end-1),xt'), wt, xt') ];
% 
% CUM_rc = CumIG(wt,update(end),mu(update(1:end-1),xt'));
% 
% IG_rc =InvGauss(wt,mu(update(1:end-1),xt'),update(end));
% 
% GRAD_RC = (grad_IG_rc .* (1-CUM_rc) - IG_rc .* grad_CUM_rc) ./ (1-CUM_rc).^2 ;

%%

HESS=zeros(length(update));

for i=1:length(update)
    h_i=zeros(length(update),1)';
    incr_i=nthroot(eps,3)*(1+abs(update(i)));
    h_i(i)=incr_i;

    incr_i=update+h_i;
    
    grad_LogIG_incr = [ dthLogIG(incr_i(end),mu(incr_i(1:end-1),xn),wn,xn) , dkLogIG(incr_i(end),mu(incr_i(1:end-1),xn),wn) ];
    grad_CUM_incr = [ dkCumIG(incr_i(end),mu(incr_i(1:end-1),xn), wn) , dthCumIG(incr_i(end),mu(incr_i(1:end-1),xn), wn, xn) ];
    CUM_incr = CumIG(wn,incr_i(end),mu(incr_i(1:end-1),xn));
    GRAD_U_incr = (grad_LogIG_incr + grad_CUM_incr ./ (1-CUM_incr)) ;
    GRAD_U_incr=sum(GRAD_U_incr .* 0.005);

    HESS(:,i)=  (GRAD_U_incr - GRAD_U) ./ incr_i;

    end


FIX = GRAD/HESS;
update = update - FIX;
step = step + 1;


