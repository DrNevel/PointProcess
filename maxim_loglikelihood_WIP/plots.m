load('EKGR.mat')
L=L;
delta=0.005;
[KSdistance,Z,T,ordered,KSoutPerc,lin,lu,ll] = ks_plot(EKGR, L, delta, 1);


figure,
plot(Mu_cif)
hold on 
plot(Mu_rc)
legend('prof uncensored','CIF Uncensored')
title(['Mu CIF uncensored vs Mu IG uncensored'])

figure,
plot(Kappa_prof)
hold on 
plot(Kappa)
legend('CIF uncensored','prof Uncensored')
title(['K CIF uncensored vs K IG uncensored'])
    
figure,
i=1;
plot(Thetap_cif(i,:))
hold on 
plot(Thetap(i+1,:))
legend('CIF uncensored','prof Uncensored')
title(['K(',  num2str(i),') CIF uncensored vs IG uncensored'])
