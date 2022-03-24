load('EKGR.mat')
L=L_cif_rc;
delta=0.005;
[KSdistance,Z,T,ordered,KSoutPerc,lin,lu,ll] = ks_plot(EKGR, L, delta, 1);


figure,
plot(Mu_cif)
hold on 
plot(Mu_cif_rc)
legend('CIF u','CIF rc')

figure,
plot(Kappa_prof)
hold on 
plot(Kappa)
legend('CIF rc','prof rc')
    
figure,
i=1;
plot(Thetap_cif(i,:))
hold on 
plot(Thetap(i+1,:))
legend('CIF rc','prof rc')

