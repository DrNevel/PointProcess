function [PP] = get_PP(EKGR, varargin)

% Copyright (C) Maximiliano Mollura and Riccardo Barbieri, 2019-2020.
% All Rights Reserved. See LICENSE.TXT for license details.
% maximiliano.mollura@polimi.it
% barbieri@neurostat.mit.edu

delta = .005; % seconds %% Compute at 200 Hz
UndSampFact = 1; % Save at 1/(UndSampFact*delta) % 20 Hz
P = 9;
W = 90;
get_pparrythmia = 0;
get_mono = 0; COV_T = []; COV_VAL = [];
get_multi = 0;
get_biva = 0;
% get_spectra = 0;
rrcorrection = [];

for i = 1:2:length(varargin)
    switch varargin{i}
        case 'delta'
            delta = varargin{i+1};
        case 'UndSampl'
            UndSampFact = varargin{i+1};
        case 'order'
            P  = varargin{i+1};
        case 'W'
            W = varargin{i+1};
        case 'get_mono'
            get_mono = varargin{i+1};
            get_biva = 0;
            get_multi = 0;
        case 'cov_t'
            COV_T = varargin{i+1};
        case 'cov_v'
            COV_VAL = varargin{i+1};
        case 'get_biva'
            get_biva = varargin{i+1};
            get_mono = 0;
            get_multi = 0;
        case 'get_multi'
            get_multi = varargin{i+1};
            get_mono = 0;
            get_biva = 0;
            %         case 'get_spectra'
            %             get_spectra = varargin{i+1};
        case 'pparrhy'
            get_pparrythmia = varargin{i+1};
        case 'rr_corr'
            rrcorrection = varargin{i+1};
        case 'Pmono'
            Pmono = varargin{i+1};
    end
end

try
    
    if get_mono
        %% MONOVARIATE PP MODEL
        if get_pparrythmia
            %             disp(['PP arrythmia...' newline])
            %             try
            %                 cov_str = {'SYST' 'SYST_VAL'};
            %                 [EKGR_new,Action,~,opt_arr] = pparrythmia_V1(EKGR(:), 'hasTheta0', 1, 'P', P,'W',W ...
            %                     , 'covariates', {[SYST(:), SYST_VAL(:)]},'toPlot',0);
            %             catch ME
            %                 disp(['PP arrythmia not working' newline])
            %                 disp(ME.message)
            %             end
        end
        
        disp(['Building Monovariate Point Process model...' newline])
        
        [Thetap,Mu,Kappa,L,opt] = pplikel_corretto(EKGR(:), 'hasTheta0', 1, 'delta', delta, 'P', P,'W',W);
               
        % KS-Distance
        [KSdistance,Z,T,ordered,KSoutPerc,lin,lu,ll] = ks_plot(EKGR, L, delta, 0);
        % ACF
        %         [acf,lags,bounds] = autocorr(T,'NumLags',60);
        [acf,bounds] = check_corr(Z,60,0);
        acf_over_thr = sum(abs(acf) > bounds((bounds>0)));
        % Variance
        if ~isempty(rrcorrection) % Correcting the Thetap we can have the correct length for all the spectral parameters
            % rrcorrection = [seconds of the event, duration of the event in seconds, replaced value in seconds]
            rrcorrection(:,1) = rrcorrection(:,1) - opt.t0; % depolarization
            cum_rrcorr = [0; cumsum(rrcorrection(1:end-1,2)-rrcorrection(1:end-1,3))];
            [Mu,Thetap,Kappa,L,opt.meanRR,opt.LogLikel] = openECG(rrcorrection,cum_rrcorr,delta,Mu,Thetap,Kappa,L,opt.meanRR,opt.LogLikel);
        end
        
        s2 = opt.meanRR.^3./Kappa;
        
        
    elseif get_biva
        %% BIVARIATE PP MODEL
        if get_pparrythmia
            %                 disp(['PP arrythmia...' newline])
            %                 try
            %                     cov_str = {'SYST' 'SYST_VAL'};
            %                     [EKGR_new,Action,~,opt_arr] = pparrythmia_V1(EKGR(:), 'hasTheta0', 1, 'P', P,'W',W ...
            %                         , 'covariates', {[SYST(:), SYST_VAL(:)]},'toPlot',0);
            %                 catch ME
            %                     disp(['PP arrythmia not working' newline])
            %                     disp(ME.message)
            %                 end
        end
        
        disp(['Building Bivariate Point Process model...' newline])
        
        if exist('Pmono','var') %% Compute also Monovariate on the Covariates (useful for Granger Causality Index)
            [Thetapc,Muc,Kappac,Lc,optc] = pplikel_cov_corr(EKGR, {[COV_T(:), COV_VAL(:)]},...
                'P', P, 'hasTheta0', 1, 'delta', delta,'W', W, 'P2',Pmono);
        else
            [Thetapc,Muc,Kappac,Lc,optc] = pplikel_cov_corr(EKGR, {[COV_T(:), COV_VAL(:)]},...
                'P', P, 'hasTheta0', 1, 'delta', delta,'W', W);
        end
        
        disp(newline)
              
        % KS-Distance
        [KSdistancec,Zc,Tc,orderedc,KSoutPercCov,linc,luc,llc] = ks_plot(EKGR, Lc, delta, 0);
        % ACF
        [acf,bounds] = check_corr(Zc,60,0);
        acf_over_thr = sum(abs(acf) > bounds((bounds>0)));
        acf_5samp_over_thr = sum(abs(acf(1:5)) > bounds((bounds>0)));
        
        % Realign in Time if closed
        if ~isempty(rrcorrection)
            rrcorrection(:,1) = rrcorrection(:,1) - optc.t0; % depolarization
            cum_rrcorr = [0; cumsum(rrcorrection(1:end-1,2)-rrcorrection(1:end-1,3))];
            [Muc,Thetapc,Kappac,Lc,optc.meanRR,~,optc.Thetap2,optc.Var2] = ...
                openECG(rrcorrection,cum_rrcorr,delta,Muc,Thetapc,Kappac,Lc,optc.meanRR,[],optc.Thetap2,optc.Var2);
        end
        % Variances
        s2c = optc.meanRR.^3./Kappac;
        cov_s2c = optc.Var2; % covariates variance
        
        
    elseif get_multi
        %% MULTIVARIATE PP MODEL
        % Coming Soon...
        
    end
    %% SAVING
    disp('Exiting PP Modeling...')
    % Features Saving
    if get_mono
        % Monovariate
        PP.Mu = Mu(:,1:UndSampFact:end); clear Mu
        PP.Thetap = Thetap(:,1:UndSampFact:end); clear Thetap
        PP.meanRR = opt.meanRR(:,1:UndSampFact:end);
        PP.Kappa = Kappa(:,1:UndSampFact:end); clear Kappa
        PP.L = L(:,1:UndSampFact:end); clear L
        PP.s2 = s2(:,1:UndSampFact:end); clear s2
        PP.delta = delta*UndSampFact;
        PP.t0 = opt.t0;
        PP.LogLikel = opt.LogLikel(:,1:UndSampFact:end); clear opt
        PP.KSdistance = KSdistance; clear KSdistance
        PP.KSoutPerc = KSoutPerc; clear KSoutPercCov
        PP.Z = Z; clear Z
        PP.T = T; clear T
        PP.acf = acf; clear acf
        PP.bounds = bounds; clear bounds
        PP.acf_over_thr = acf_over_thr; clear acf_over_thr
        
    end
    if get_biva
        % Covariate
        PP.Mu = Muc(:,1:UndSampFact:end);  clear Muc
        PP.Thetap = Thetapc(:,1:UndSampFact:end); clear Thetapc
        PP.meanRR = optc.meanRR(:,1:UndSampFact:end);
        PP.Kappa = Kappac(:,1:UndSampFact:end); clear Kappac
        PP.L = Lc(:,1:UndSampFact:end); clear Lc
        PP.s2 = s2c(:,1:UndSampFact:end); clear s2c
        PP.cov_s2 = cov_s2c(:,1:UndSampFact:end); clear cov_s2c
        PP.delta = delta*UndSampFact;
        PP.t0 = optc.t0;
        if isfield(optc,'Thetap_cov')
            PP.Thetap_cov = optc.Thetap_cov(:,1:UndSampFact:end,:);
            PP.Var_cov = optc.Var_Cov(:,1:UndSampFact:end);
        end
        PP.Thetap2 = optc.Thetap2(:,1:UndSampFact:end,:); clear optc
        PP.KSdistance = KSdistancec; clear KSdistancec 
        PP.KSoutPerc = KSoutPercCov; clear KSoutPercCov
        PP.Z = Zc; clear Zc
        PP.T = Tc; clear Tc
        PP.acf = acf; clear acf
        PP.bounds = bounds; clear bounds
        PP.acf_over_thr = acf_over_thr; clear acf_over_thr
    end
    
%     if get_multi
%         PP.Mu = Muc(:,1:UndSampFact:end);  clear Muc
%         PP.Thetap = Thetapc(:,1:UndSampFact:end); clear Thetapc
%         PP.meanRR = optc.meanRR(:,1:UndSampFact:end);
%         PP.Kappa = Kappac(:,1:UndSampFact:end); clear Kappac
%         PP.L = Lc(:,1:UndSampFact:end); clear Lc
%         PP.tau = tau;
%         PP.s2 = s2c(:,1:UndSampFact:end); clear s2c
%         PP.cov_s2 = cov_s2c(:,1:UndSampFact:end); clear cov_s2c
%         PP.delta = delta*UndSampFact;
%         PP.t0 = optc.t0;
%         PP.covariates_ar = optc.Thetap2(:,:,1:UndSampFact:end); clear optc
%         PP.KSdistance = KSdistancec; clear KSdistancec
%         PP.Z = Zc; clear Zc
%         PP.T = Tc; clear Tc
%         PP.acf = acf; clear acf
%         PP.bounds = bounds; clear bounds
%         PP.acf_over_thr = acf_over_thr; clear acf_over_thr
% 
%     end
catch ME
    disp('Failed')
    disp(ME.message)
    PP = deal([]);
    PP.error = ME;
end

end