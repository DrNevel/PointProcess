function GRAD = eval_grad(CumIG, dkLogIG, dkCumIG, dthLogIG, dthCumIG, mu, update, wn, xn,eta)
    kLogIG= dkLogIG(update(end),mu(update(1:end-1),xn),wn);
    thLogIG= dthLogIG(update(end),mu(update(1:end-1),xn),wn,xn);
    kCumIG=dkCumIG(update(end),mu(update(1:end-1),xn), wn);
    thCumIG=dthCumIG(update(end),mu(update(1:end-1),xn), wn, xn);
    CIG=CumIG(wn,update(end),mu(update(1:end-1),xn));
    
    

    GRAD = sum(eta.*[ thLogIG , kLogIG ]) +  sum(eta.*[ thCumIG , kCumIG ]) ./ sum(eta.*(1-CIG));
end
