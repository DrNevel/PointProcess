function GRAD = eval_grad(CumIG, dkLogIG, dkCumIG, dthLogIG, dthCumIG, mu, update, wn, xn)
    kLogIG= dkLogIG(update(end),mu(update(1:end-1),xn),wn);
    thLogIG= dthLogIG(update(end),mu(update(1:end-1),xn),wn,xn);
    kCumIG=dkCumIG(update(end),mu(update(1:end-1),xn), wn);
    thCumIG=dthCumIG(update(end),mu(update(1:end-1),xn), wn, xn);
    CIG=CumIG(wn,update(end),mu(update(1:end-1),xn));
    
    

    GRAD = sum([ thLogIG , kLogIG ]) +  sum([ thCumIG , kCumIG ]) ./ sum(1-CIG);
end
