function GRAD_rc = eval_grad_rc(CumIG, dkIG, dkCumIG, dthIG, dthCumIG, InvGauss, mu, update, wt, xt)
    
%     WT=wt;

    if isscalar(wt)
        WT=[0.005:0.005:wt];
    else
        WT=wt;
    end

    kIG= dkIG(update(end),mu(update(1:end-1),xt'),WT);
    thIG= dthIG(update(end),mu(update(1:end-1),xt'),WT,xt);
    kCumIG=dkCumIG(update(end),mu(update(1:end-1),xt'), WT);
    thCumIG=dthCumIG(update(end),mu(update(1:end-1),xt'), WT, xt);
    IG = InvGauss(WT, mu(update(1:end-1),xt'), update(end));
    CIG=CumIG(WT,update(end),mu(update(1:end-1),xt'));
    
   
% Gradient with scalar wt    
%     GRAD_rc = ( ( ([ thIG' , kIG ]).*(1-CIG) +  (IG).*([ thCumIG' , kCumIG ]) ) ./ (0.99*(1-CIG).^2) ) *wt  ;

% Gradient with vector wt
    GRAD_rc = ( ( ([thIG ; kIG]).*(1-CIG) + (IG).*([thCumIG ; kCumIG]) ) ./ (0.99*(1-CIG).^2) ) ;
end
