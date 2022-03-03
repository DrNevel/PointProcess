function [AIC_v,FPE_v,PEWT_v] = AR_assess(sig,orders)

k = 1;
for i = orders
    
    modl=ar(sig,i,'yw');
    AIC_v(k) = aic(modl);
    FPE_v(k) = fpe(modl);
    [e,v]=resid(sig(:),modl);
    PEWT_v(k)=adtest(e.OutputData);
    k = k+1;
end
