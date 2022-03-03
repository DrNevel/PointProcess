function [Mu,Thetap,Theta0,Kappa,L,meanRR,LogLikel,covariates_ar,covariates_s2,Mu_Cov,Thetap_mono_cov,Var_mono_cov,Mean_cov] = openECG(rrcorrection,cum_rrcorr,delta,Mu,Thetap,Theta0,...
    Kappa,L,meanRR,LogLikel,covariates_ar,covariates_s2,Mu_Cov, Thetap_mono_cov,Var_mono_cov,Mean_cov, covariates_ar_t,covariates_s2_t)

% Copyright (C) Maximiliano Mollura and Riccardo Barbieri, 2019-2020.
% All Rights Reserved. See LICENSE.TXT for license details.
% maximiliano.mollura@polimi.it
% riccardo.barbieri@polimi.it


% warning('Check if the opening is working properly!')

    if ~exist('Mu','var')
        Mu = [];
    end
    if ~exist('Thetap','var')
        Thetap = [];
    end
    if ~exist('Theta0','var')
        Theta0 = [];
    end
    if ~exist('Kappa','var')
        Kappa = [];
    end
    if ~exist('L','var')
        L = [];
    end
    if ~exist('meanRR','var')
        meanRR = [];
    end
    if ~exist('LogLikel','var')
        LogLikel = [];
    end
    if ~exist('covariates_ar','var')
        covariates_ar = [];
    end
    if ~exist('covariates_s2','var')
        covariates_s2 = [];
    end
    if ~exist('Mu_Cov','var')
        Mu_Cov = [];
    end
    if ~exist('covariates_ar_t','var')
        covariates_ar_t = [];
    end
    if ~exist('covariates_s2_t','var')
        covariates_s2_t = [];
    end
    if ~exist('Thetap_mono_cov','var')
        Thetap_mono_cov = [];
    end 
    if ~exist('Var_mono_cov','var')
        Var_mono_cov = [];
    end
    if ~exist('Mean_cov','var')
        Mean_cov = [];
    end
    
for i = 1:size(rrcorrection,1)
    if ~isempty(Mu)
        tmp = [Mu(:,1:round((rrcorrection(i,1)+cum_rrcorr(i))/delta)),NaN(size(Mu,1),round((rrcorrection(i,2)-rrcorrection(i,3))/delta)),Mu(:,round((rrcorrection(i,1)+cum_rrcorr(i))/delta)+1:end)];
        Mu = tmp;
    end
    if ~isempty(Thetap)
        tmp = [Thetap(:,1:round((rrcorrection(i,1)+cum_rrcorr(i))/delta)),NaN(size(Thetap,1),round((rrcorrection(i,2)-rrcorrection(i,3))/delta)),Thetap(:,round((rrcorrection(i,1)+cum_rrcorr(i))/delta)+1:end)];
        Thetap = tmp;
    end
    if ~isempty(Theta0)
        tmp = [Theta0(:,1:round((rrcorrection(i,1)+cum_rrcorr(i))/delta)),NaN(size(Theta0,1),round((rrcorrection(i,2)-rrcorrection(i,3))/delta)),Theta0(:,round((rrcorrection(i,1)+cum_rrcorr(i))/delta)+1:end)];
        Theta0 = tmp;
    end    
    if ~isempty(Kappa)
        tmp = [Kappa(:,1:round((rrcorrection(i,1)+cum_rrcorr(i))/delta)),NaN(size(Kappa,1),round((rrcorrection(i,2)-rrcorrection(i,3))/delta)),Kappa(:,round((rrcorrection(i,1)+cum_rrcorr(i))/delta)+1:end)];
        Kappa = tmp;
    end
    if ~isempty(L)
        tmp = [L(:,1:round((rrcorrection(i,1)+cum_rrcorr(i))/delta)),NaN(size(L,1),round((rrcorrection(i,2)-rrcorrection(i,3))/delta)),L(:,round((rrcorrection(i,1)+cum_rrcorr(i))/delta)+1:end)];
        L = tmp;
    end
    if ~isempty(meanRR)
        tmp = [meanRR(:,1:round((rrcorrection(i,1)+cum_rrcorr(i))/delta)),NaN(size(meanRR,1),round((rrcorrection(i,2)-rrcorrection(i,3))/delta)),meanRR(:,round((rrcorrection(i,1)+cum_rrcorr(i))/delta)+1:end)];
        meanRR = tmp;
    end
    if ~isempty(LogLikel)
        tmp = [LogLikel(:,1:round((rrcorrection(i,1)+cum_rrcorr(i))/delta)),NaN(size(LogLikel,1),round((rrcorrection(i,2)-rrcorrection(i,3))/delta)),LogLikel(:,round((rrcorrection(i,1)+cum_rrcorr(i))/delta)+1:end)];
        LogLikel = tmp;
    end
    if ~isempty(covariates_ar)
        tmp = [covariates_ar(:,1:round((rrcorrection(i,1)+cum_rrcorr(i))/delta),:),NaN(size(covariates_ar,1),round((rrcorrection(i,2)-rrcorrection(i,3))/delta),size(covariates_ar,3)),covariates_ar(:,round((rrcorrection(i,1)+cum_rrcorr(i))/delta)+1:end,:)];
        covariates_ar = tmp;
    end
    if ~isempty(covariates_s2)
        tmp = [covariates_s2(:,1:round((rrcorrection(i,1)+cum_rrcorr(i))/delta)),NaN(size(covariates_s2,1),round((rrcorrection(i,2)-rrcorrection(i,3))/delta)),covariates_s2(:,round((rrcorrection(i,1)+cum_rrcorr(i))/delta)+1:end)];
        covariates_s2 = tmp;
    end
    if ~isempty(Mu_Cov)
        tmp = [Mu_Cov(:,1:round((rrcorrection(i,1)+cum_rrcorr(i))/delta)),NaN(size(Mu_Cov,1),round((rrcorrection(i,2)-rrcorrection(i,3))/delta)),Mu_Cov(:,round((rrcorrection(i,1)+cum_rrcorr(i))/delta)+1:end)];
        Mu_Cov = tmp;
    end
    if ~isempty(Thetap_mono_cov)
        tmp = [Thetap_mono_cov(:,1:round((rrcorrection(i,1)+cum_rrcorr(i))/delta),:),NaN(size(Thetap_mono_cov,1),round((rrcorrection(i,2)-rrcorrection(i,3))/delta),size(Thetap_mono_cov,3)),Thetap_mono_cov(:,round((rrcorrection(i,1)+cum_rrcorr(i))/delta)+1:end,:)];
        Thetap_mono_cov = tmp;
    end
    if ~isempty(Var_mono_cov)
        tmp = [Var_mono_cov(:,1:round((rrcorrection(i,1)+cum_rrcorr(i))/delta)),NaN(size(Var_mono_cov,1),round((rrcorrection(i,2)-rrcorrection(i,3))/delta),size(Var_mono_cov,1)),Var_mono_cov(:,round((rrcorrection(i,1)+cum_rrcorr(i))/delta)+1:end)];
        Var_mono_cov = tmp;
    end
    if ~isempty(Mean_cov)
        tmp = [Mean_cov(:,1:round((rrcorrection(i,1)+cum_rrcorr(i))/delta)),NaN(size(Mean_cov,1),round((rrcorrection(i,2)-rrcorrection(i,3))/delta),size(Mean_cov,1)),Mean_cov(:,round((rrcorrection(i,1)+cum_rrcorr(i))/delta)+1:end)];
        Mean_cov = tmp;
    end
end


end