GRAD_rc_tot=0;

for i=1:length(wn)
    GRAD_temp = eval_grad_rc(CumIG, dkIG, dkCumIG, dthIG, dthCumIG, InvGauss, mu, update, wn(i), xn);

    GRAD_rc_tot = GRAD_rc_tot + GRAD_temp;
end