function [PP_new, PP_U, PP] = get_PP_edit(EKGR, varargin)

% Copyright (C) Maximiliano Mollura and Riccardo Barbieri, 2019-2020.
% All Rights Reserved. See LICENSE.TXT for license details.
% maximiliano.mollura@polimi.it
% riccardo.barbieri@polimi.it

% last edit: 19/10/2020 - Added covariates to Mu
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
                EKGR_or=EKGR;
            try
                EKGR=pparrythmia(EKGR,'hasTheta0', 1, 'delta', delta, 'P', P,'W',W);
            catch ME
                disp(['PP arrythmia not working' newline])
                disp(ME.message)
            end
        end
        
        disp(['Building Monovariate Point Process model...' newline])
        
        %[Thetap,Mu,Kappa,L,opt] = pplikel_corretto(EKGR(:), 'hasTheta0', 1, 'delta', delta, 'P', P,'W',W);

        
        % RC new (non funziona)
        [Thetap_rc,Mu_rc,Kappa_rc,L_rc,opt_rc] = pplikel_new_rc(EKGR(:), 'hasTheta0', 1, 'delta', delta, 'P', P,'W',W);

        % Uncensored new (nostro)
        [Thetap_new,Mu_new,Kappa_new,LogLikel_new,opt_new] = pplikel_new(EKGR(:), 'hasTheta0', 1, 'delta', delta, 'P', P,'W',W);
        
        % RC prof
        [Thetap,Mu,Kappa,L,opt] = pplikel_corretto(EKGR(:), 'hasTheta0', 1, 'delta', delta, 'P', P,'W',W);

        % Uncensored prof (modificata da noi)
        [Thetap_U,Mu_U,Kappa_U,opt_U] = pplikel_corretto_uncensored(EKGR(:), 'hasTheta0', 1, 'delta', delta, 'P', P,'W',W);
                
          
        %%%%% PLOT MU
%         figure        
%         plot(Mu_U);
%         hold on
%         plot(Mu_new);
%         title('Superimposed Mu')
%         legend('Mu Uncensored','Mu New Estimate')
%         
%         %%%% PLOT SHAPE FACTOR K
%         figure
%         plot(Kappa_U);
%         hold on
%         plot(Kappa_new);
%         title('Superimposed Shape Factor')
%         legend('K uncensored', 'Kappa new estimate')
        
        %%%% PLOT THETA
%         for i=1:size(Thetap_new,1)
%             figure
%             plot(Thetap_U(i,:));
%             hold on
%             plot(Thetap_new(i,:));
%             legend('Theta Uncensored','Theta New Estimate')
%             title(['Theta(',  num2str(i),') uncensored vs new estimate'])
%         end
        
               
        
        % KS-Distance (compute only in case we have L) --> ci manca L per
        % la uncensored old e la uncensored new
       [KSdistance,Z,T,ordered,KSoutPerc,lin,lu,ll] = ks_plot(EKGR, L, delta, 0);

       % ACF
       %         [acf,lags,bounds] = autocorr(T,'NumLags',60);
       [acf,bounds] = check_corr(Z,60,0);
       acf_over_thr = sum(abs(acf) > bounds((bounds>0)));
                
        
        % Variances
        if ~isempty(rrcorrection) % Correcting the Thetap we can have the correct length for all the spectral parameters
            % rrcorrection = [seconds of the event, duration of the event in seconds, replaced value in seconds]
            rrcorrection(:,1) = rrcorrection(:,1) - opt.t0; % depolarization
            cum_rrcorr = [0; cumsum(rrcorrection(1:end-1,2)-rrcorrection(1:end-1,3))];
            % Correction NEW
            [Mu_new,Thetap_new,opt_new.Theta0,Kappa_new,opt_new.meanRR,LogLikel_new] = openECG(rrcorrection,cum_rrcorr,delta,Mu_new,Thetap_new,opt_new.Theta0,Kappa_new,opt_new.meanRR,LogLikel_new);
            % Correction OLD UNCENSORED
            [Mu_U,Thetap_U,opt_U.Theta0,Kappa_U,opt_U.meanRR] = openECG(rrcorrection,cum_rrcorr,delta,Mu_U,Thetap_U,opt_U.Theta0,Kappa_U,opt_U.meanRR);
            % Correction CENSORED
            [Mu,Thetap,opt.Theta0,Kappa,L,opt.meanRR,opt.LogLikel] = openECG(rrcorrection,cum_rrcorr,delta,Mu,Thetap,opt.Theta0,Kappa,L,opt.meanRR,opt.LogLikel);
        end
                  
        s2 = opt.meanRR.^3./Kappa;
        s2_U = opt_U.meanRR.^3./Kappa_U;
        s2_new = opt_new.meanRR.^3./Kappa_new;
        
        %%%% PLOT SD
%         figure
%         plot(s2_U);
%         hold on
%         plot(s2_new);
%         title('Superimposed Sstandard Deviatiom')
%         legend('SD uncensored', 'SD new estimate')
        
        
    elseif get_biva
        %% BIVARIATE PP MODEL
        % funziona per il biva.....?
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
            [Muc,Thetapc,optc.Theta0,Kappac,Lc,optc.meanRR,~,optc.Thetap2,optc.Var2,optc.Cov_Gaus,optc.meanCOV] = ...
                openECG(rrcorrection,cum_rrcorr,delta,Muc,Thetapc,optc.Theta0,Kappac,Lc,optc.meanRR,[],optc.Thetap2,optc.Var2,optc.Cov_Gaus,optc.meanCOV);
        end
        % Variances
        s2c = optc.meanRR.^3./Kappac;
        cov_s2c = optc.Var2; % covariates variance
        
        
    elseif get_multi
        %% MULTIVARIATE PP MODEL
        [Thetapc,Muc,Kappac,Lc,optc] = pplikel_cov_corr_multi(EKGR, {[COV_T(1,:);COV_VAL(1,:)],[COV_T(2,:);COV_VAL(2,:)]},...
            'P', P, 'hasTheta0', 1, 'delta', delta,'W', W);
        
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
            [Muc,Thetapc,optc.Theta0,Kappac,Lc,optc.meanRR,~,optc.Thetap2,optc.Var2,optc.Cov_Gaus,optc.meanCOV] = ...
                openECG(rrcorrection,cum_rrcorr,delta,Muc,Thetapc,optc.Theta0,Kappac,Lc,optc.meanRR,[],optc.Thetap2,optc.Var2,optc.Cov_Gaus,optc.meanCOV);
        end
        % Variances
        s2c = optc.meanRR.^3./Kappac;
        cov_s2c = optc.Var2; % covariates variance

        % Coming Soon...
        
    end
    %% SAVING
    disp('Exiting PP Modeling...')
    % Features Saving
    if get_mono
        % Monovariate NEW
        PP_new.Mu = Mu_new(:,1:UndSampFact:end); clear Mu_new 
        PP_new.Thetap = Thetap_new(:,1:UndSampFact:end); clear Thetap_new
        if opt_new.hasTheta0,PP_new.Theta0 = opt_new.Theta0; end
        PP_new.meanRR = opt_new.meanRR(:,1:UndSampFact:end);
        PP_new.Kappa = Kappa_new(:,1:UndSampFact:end); clear Kappa_new
        %PP_new.L = L(:,1:UndSampFact:end); clear L_new
        PP_new.s2 = s2_new(:,1:UndSampFact:end); clear s2_new
        PP_new.delta = delta*UndSampFact;
        PP_new.t0 = opt_new.t0;
        PP_new.LogLikel = LogLikel_new(:,1:UndSampFact:end);
        PP_new.L=opt_new.L;
        clear opt_new
        %PP_new.KSdistance = KSdistance; clear KSdistance
        %PP_new.KSoutPerc = KSoutPerc; clear KSoutPercCov
        %PP_new.Z = Z; clear Z
        %PP_new.T = T; clear T
        %PP_new.acf = acf; clear acf
        %PP_new.bounds = bounds; clear bounds
        %PP_new.acf_over_thr = acf_over_thr; clear acf_over_thr
        
        %%%% MONOVARIATE OLD UNCENSORED
        PP_U.Mu = Mu_U(:,1:UndSampFact:end); clear Mu_U
        PP_U.Thetap = Thetap_U(:,1:UndSampFact:end); clear Thetap_U
        if opt_U.hasTheta0,PP_U.Theta0 = opt_U.Theta0; end
        PP_U.meanRR = opt_U.meanRR(:,1:UndSampFact:end);
        PP_U.Kappa = Kappa_U(:,1:UndSampFact:end); clear Kappa_U
        %PP_U.L = L(:,1:UndSampFact:end); clear L_U
        PP_U.s2 = s2_U(:,1:UndSampFact:end); clear s2_U
        PP_U.delta = delta*UndSampFact;
        PP_U.t0 = opt_U.t0;
        %PP_U.LogLikel = LogLikel_U(:,1:UndSampFact:end); clear opt_U
        %PP_new.KSdistance = KSdistance; clear KSdistance
        %PP_new.KSoutPerc = KSoutPerc; clear KSoutPercCov
        %PP_new.Z = Z; clear Z
        %PP_new.T = T; clear T
        %PP_new.acf = acf; clear acf
        %PP_new.bounds = bounds; clear bounds
        %PP_new.acf_over_thr = acf_over_thr; clear acf_over_thr
        
        
        % MONOVARIATE NEW RIGHT CENSORED
%         PP_rc.Mu = Mu_rc(:,1:UndSampFact:end); clear Mu 
%         PP_rc.Thetap = Thetap_rc(:,1:UndSampFact:end); clear Thetap
%         if opt_rc.hasTheta0,PP_rc.Theta0 = opt_rc.Theta0; end
%         PP_rc.meanRR = opt_rc.meanRR_rc(:,1:UndSampFact:end);
%         PP_rc.Kappa = Kappa_rc(:,1:UndSampFact:end); clear Kappa
%         PP_rc.L = L_rc(:,1:UndSampFact:end); clear L
%         PP_rc.s2 = s2_rc(:,1:UndSampFact:end); clear s2
%         PP_rc.delta = delta*UndSampFact;
%         PP_rc.t0 = opt_rc.t0;
%         PP_rc.LogLikel = opt_rc.LogLikel_rc(:,1:UndSampFact:end); clear opt
        
        %PP_rc.KSdistance = KSdistance; clear KSdistance
        %PP_rc.KSoutPerc = KSoutPerc; clear KSoutPercCov
        %PP_rc.Z = Z; clear Z
        %PP_rc.T = T; clear T
        %PP_rc.acf = acf; clear acf
        %PP_rc.bounds = bounds; clear bounds
        %PP_rc.acf_over_thr = acf_over_thr; clear acf_over_thr
        
        
        % MONOVARIATE CENSORED (prof)
        PP.Mu = Mu(:,1:UndSampFact:end); clear Mu 
        PP.Thetap = Thetap(:,1:UndSampFact:end); clear Thetap
        if opt.hasTheta0,PP.Theta0 = opt.Theta0; end
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
        PP.Mu = [Muc(:,1:UndSampFact:end); optc.Cov_Gaus(:,1:UndSampFact:end)];  clear Muc 
        PP.Thetap = Thetapc(:,1:UndSampFact:end); clear Thetapc
        PP.meanRR = optc.meanRR(:,1:UndSampFact:end);
        PP.meanCOV = optc.meanCOV(:,1:UndSampFact:end);
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
    
    if get_multi
        % Covariate
        PP.Mu = [Muc(:,1:UndSampFact:end); optc.Cov_Gaus(:,1:UndSampFact:end)];  clear Muc 
        PP.Thetap = Thetapc(:,1:UndSampFact:end); clear Thetapc
        PP.meanRR = optc.meanRR(:,1:UndSampFact:end);
        PP.meanCOV = optc.meanCOV(:,1:UndSampFact:end);
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
catch ME
    disp('Failed')
    disp(ME.message)
    PP = deal([]);
    PP.error = ME;
end

end