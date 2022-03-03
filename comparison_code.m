%% load paramters
clear; 

stored=load('param_RC_NEW_2');
PP_rc_new.Thetap=stored.Thetap;
PP_rc_new.Mu=stored.Mu;
PP_rc_new.K=stored.Kappa;
PP_rc_new.meanRR=stored.meanRR;
PP_rc_new.Loglikel=stored.Loglikel;
PP_rc_new.opt=stored.opt;
clear 'stored';

stored=load('param_RC_OLD');
PP_rc_old.Thetap=stored.Thetap;
PP_rc_old.Mu=stored.Mu;
PP_rc_old.K=stored.Kappa;
PP_rc_old.L=stored.L;
PP_rc_old.Loglikel=stored.opt.LogLikel;
PP_rc_old.opt=stored.opt;
clear 'stored';

stored=load('param_U_NEW');
PP_u_new.Thetap=stored.Thetap_new;
PP_u_new.Mu=stored.Mu_new;
PP_u_new.K=stored.Kappa_new;
PP_u_new.opt=stored.opt_new;
PP_u_new.Loglikel=stored.LogLikel_new;
clear 'stored';


stored=load('parameters_from_estim');
estim.L=stored.opt.L_estim;
estim.LogLikel=stored.opt.LL_estim;
estim.IG=stored.opt.IG_estim;
estim.CIG=stored.opt.CIG_estim;
clear 'stored';

%% load GRAD, HESS, LL, Lambda
stored=load('workspace_comparison');
GRAD.RC.it1=stored.GRAD_it1_rc;
GRAD.U.itlast=stored.GRAD_itlast_u;
GRAD.U.it1=stored.GRAD_it1_u;

HESS.RC.it1=stored.HESS_it1_rc;
HESS.U.itlast=stored.HESS_itlast_u;
HESS.U.it1=stored.HESS_it1_u;

clear 'stored';

%% Comparisons
%% LOG LIKELIHOOD 
%%%%% old RC vs new RC
figure
plot(PP_rc_new.Loglikel)
hold on
plot(PP_rc_old.opt.LogLikel)
legend('new LL', 'old LL')

% old RC vs estimate
figure
title('LogLikel RC old vs estim')
plot(estim.LogLikel)
hold on
plot(PP_rc_old.opt.LogLikel)
legend('estim LL', 'old LL')

%%%%% new RC vs new U
figure
title('LogLikel RC old vs RC new')
plot(PP_rc_new.Loglikel)
hold on
plot(PP_u_new.Loglikel) 
legend('new LL', 'new LL')

%% HAZARD
%%%%% old RC vs new RC
figure
title('Lambda RC old vs RC new')
plot(PP_rc_new.opt.L)
hold on
plot(PP_rc_old.L)
legend('new Lambda', 'old Lambda')

%%%%% old RC vs estimate
figure
title('Lambda RC old vs estimate')
plot(estim.L)
hold on
plot(PP_rc_old.L)
legend('estim Lambda', 'old Lambda')

%% RC parameters
%%%%% old RC vs new RC
i=9;
figure
title('Theta RC old vs new')
plot(PP_rc_new.Thetap(i,:))
hold on
plot(PP_rc_old.Thetap(i,:))
legend('new','old')

figure
title('Kappa RC old vs new')
plot(PP_rc_new.K)
hold on
plot(PP_rc_old.K)
legend('new','old')

figure
title('Mu RC old vs new')
plot(PP_rc_new.Mu)
hold on
plot(PP_rc_old.Mu)
legend('new','old')

%%%% new RC vs new U
i=9;
figure
title('Theta RC vs U')
plot(PP_rc_new.Thetap(i,:))
hold on
plot(PP_u_new.Thetap(i,:))
legend('new RC','new U')

figure
title('Kappa RC vs U')
plot(PP_rc_new.K)
hold on
plot(PP_u_new.K)
legend('new RC','new U')

figure
title('Mu RC vs U')
plot(PP_rc_new.Mu)
hold on
plot(PP_u_new.Mu)
legend('new RC','new U')
