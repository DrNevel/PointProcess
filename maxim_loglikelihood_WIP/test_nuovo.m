GRAD_rc_tot=[];

for i=1:length(wn)
    GRAD_temp = eval_grad_rc(CumIG, dkIG, dkCumIG, dthIG, dthCumIG, InvGauss, mu, update, wn(i), xt, eta);
    GRAD_rc_tot = [ GRAD_rc_tot , GRAD_temp] ;
       
end

GRAD_rc_tot = sum(GRAD_rc_tot.*eta',2);