function GRAD_rc = eval_grad_rc(CumIG, dkIG, dkCumIG, dthIG, dthCumIG, InvGauss, mu, update, wt, xt)
    kIG= dkIG(update(end),mu(update(1:end-1),xt'),wt);
    thIG= dthIG(update(end),mu(update(1:end-1),xt'),wt,xt);
    kCumIG=dkCumIG(update(end),mu(update(1:end-1),xt'), wt);
    thCumIG=dthCumIG(update(end),mu(update(1:end-1),xt'), wt, xt);
    IG = InvGauss(wt, mu(update(1:end-1),xt'), update(end));
    CIG=CumIG(wt,update(end),mu(update(1:end-1),xt'));
    
   
    GRAD_rc = ( ([ thIG' , kIG ]).*(1-CIG) -  (IG).*([ thCumIG' , kCumIG ]) ) ./ ((1-CIG).^2);
end

% [sum(thLogIG) , sum(kLogIG)] + [ sum(thCumIG) , sum(kCumIG) ] / sum(1-CIG)
% [sum(thCumIG./(1-CIG)) , sum(kCumIG./(1-CIG))]

% sum([ thLogIG , kLogIG ] + ( [ thCumIG , kCumIG ] ./ (1-CIG) ))