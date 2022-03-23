function GRAD_rc = eval_grad_rc(CumIG, dkIG, dkCumIG, dthIG, dthCumIG, InvGauss, mu, update, wt, xt)
    
%     WT=wt;

    WT=[0.005:0.005:wt];

    kIG= dkIG(update(end),mu(update(1:end-1),xt'),WT);
    thIG= dthIG(update(end),mu(update(1:end-1),xt'),WT,xt);
    kCumIG=dkCumIG(update(end),mu(update(1:end-1),xt'), WT);
    thCumIG=dthCumIG(update(end),mu(update(1:end-1),xt'), WT, xt);
    IG = InvGauss(WT, mu(update(1:end-1),xt'), update(end));
    CIG=CumIG(WT,update(end),mu(update(1:end-1),xt'));
    
   
    
%     GRAD_rc = ( ( ([ thIG' , kIG ]).*(1-CIG) +  (IG).*([ thCumIG' , kCumIG ]) ) ./ (0.99*(1-CIG).^2) ) *wt  ;
    GRAD_rc = ( ( ([thIG ; kIG]).*(1-CIG) + (IG).*([thCumIG ; kCumIG]) ) ./ (0.99*(1-CIG).^2) ) ;
    GRAD_rc = sum(GRAD_rc,2)*0.005;
end

% [sum(thLogIG) , sum(kLogIG)] + [ sum(thCumIG) , sum(kCumIG) ] / sum(1-CIG)
% [sum(thCumIG./(1-CIG)) , sum(kCumIG./(1-CIG))]

% sum([ thLogIG , kLogIG ] + ( [ thCumIG , kCumIG ] ./ (1-CIG) ))