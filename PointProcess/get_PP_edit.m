function [PP_prof, PP_cif] = get_PP_edit(EKGR, varargin)

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

% try
    
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
        
        % RC prof
        [Thetap_prof,Mu_prof,Kappa_prof,L_prof,opt_prof] = pplikel_corretto(EKGR(:), 'hasTheta0', 1, 'delta', delta, 'P', P,'W',W);
               
        % CIF new (uncensored funzion, rc non funziona)
        [Thetap_cif,Mu_cif,Kappa_cif,L_cif,opt_cif] = pplikel_cif(EKGR(:), 'hasTheta0', 1, 'delta', delta, 'P', P,'W',W);

        % Uncensored (nostro) ---> ???
        [Thetap_invg,Mu_invg,Kappa_invg,L_invg,opt_invg] = pplikel_invgauss(EKGR(:), 'hasTheta0', 1, 'delta', delta, 'P', P,'W',W);
%         
%         % Uncensored prof (modificata da noi)
%         [Thetap_Uprof,Mu_Uprof,Kappa_Uprof,L_Uprof,opt_Uprof] = pplikel_corretto_uncensored(EKGR(:), 'hasTheta0', 1, 'delta', delta, 'P', P,'W',W);
%         
        
        % KS-Distance (compute only in case we have L) 
       [KSdistance_prof,Z_prof,T_prof,ordered_prof,KSoutPerc_prof,lin_prof,lu_prof,ll_prof] = ks_plot(EKGR, L_prof, delta, 0);
       [KSdistance_cif,Z_cif,T_cif,ordered_cif,KSoutPerc_cif,lin_cif,lu_cif,ll_cif] = ks_plot(EKGR, L_cif, delta, 0);

       % ACF
       %         [acf,lags,bounds] = autocorr(T,'NumLags',60);
       [acf_prof,bounds_prof] = check_corr(Z_prof,60,0);
       acf_over_thr_prof = sum(abs(acf_prof) > bounds_prof((bounds_prof>0)));
       
       [acf_cif,bounds_cif] = check_corr(Z_cif,60,0);
       acf_over_thr_cif = sum(abs(acf_cif) > bounds_cif((bounds_cif>0)));
                
        
        % Variances
        if ~isempty(rrcorrection) % Correcting the Thetap we can have the correct length for all the spectral parameters
            % rrcorrection = [seconds of the event, duration of the event in seconds, replaced value in seconds]
            rrcorrection(:,1) = rrcorrection(:,1) - opt.t0; % depolarization
            cum_rrcorr = [0; cumsum(rrcorrection(1:end-1,2)-rrcorrection(1:end-1,3))];
            % Correction NEW
            [Mu_prof,Thetap_prof,opt_prof.Theta0,Kappa_prof,L_prof,opt_prof.meanRR,opt_prof.LogLikel] = openECG(rrcorrection,cum_rrcorr,delta,Mu_new,Thetap_new,opt_new.Theta0,Kappa_new,opt_new.meanRR,LogLikel_new);
            % Correction CENSORED
            [Mu_cif,Thetap_cif,opt.Theta0,Kappa_cif,L_cif,opt_cif.meanRR,opt_cif.LogLikel] = openECG(rrcorrection,cum_rrcorr,delta,Mu,Thetap,opt.Theta0,Kappa,L,opt.meanRR,opt.LogLikel);
        end
                  
        s2_prof = opt.meanRR.^3./Kappa_prof;
        s2_cif = opt_U.meanRR.^3./Kappa_cif;
        
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
        % Monovariate Prof RC
        PP_prof.Mu = Mu_prof(:,1:UndSampFact:end); clear Mu_new 
        PP_prof.Thetap = Thetap_prof(:,1:UndSampFact:end); clear Thetap_new
        if opt_new.hasTheta0,PP_prof.Theta0 = opt_new.Theta0; end
        PP_prof.meanRR = opt_prof.meanRR(:,1:UndSampFact:end);
        PP_prof.Kappa = Kappa_prof(:,1:UndSampFact:end); clear Kappa_new
        %PP_new.L = L(:,1:UndSampFact:end); clear L_new
        PP_prof.s2 = s2_prof(:,1:UndSampFact:end); clear s2_new
        PP_prof.delta = delta*UndSampFact;
        PP_prof.t0 = opt_prof.t0;
        PP_prof.LogLikel = opt_prof.LogLikel(:,1:UndSampFact:end);
        PP_prof.L=L_prof;
        clear opt_new
        %PP_new.KSdistance = KSdistance; clear KSdistance
        %PP_new.KSoutPerc = KSoutPerc; clear KSoutPercCov
        %PP_new.Z = Z; clear Z
        %PP_new.T = T; clear T
        %PP_new.acf = acf; clear acf
        %PP_new.bounds = bounds; clear bounds
        %PP_new.acf_over_thr = acf_over_thr; clear acf_over_thr
        
        %%%% MONOVARIATE CIF RC
        PP_cif.Mu = Mu_cif(:,1:UndSampFact:end); clear Mu_cif
        PP_cif.Thetap = Thetap_cif(:,1:UndSampFact:end); clear Thetap_cif
        if opt_cif.hasTheta0,PP_cif.Theta0 = opt_cif.Theta0; end
        PP_cif.meanRR = opt_cif.meanRR(:,1:UndSampFact:end);
        PP_cif.Kappa = Kappa_cif(:,1:UndSampFact:end); clear Kappa_cif
        PP_cif.L = L_cif(:,1:UndSampFact:end); clear L_cif
        PP_cif.s2 = s2_cif(:,1:UndSampFact:end); clear s2_cif
        PP_cif.delta = delta*UndSampFact;
        PP_cif.t0 = opt_cif.t0;
        PP_cif.LogLikel = opt_cif.LogLikel(:,1:UndSampFact:end); clear opt_U
        %PP_new.KSdistance = KSdistance; clear KSdistance
        %PP_new.KSoutPerc = KSoutPerc; clear KSoutPercCov
        %PP_new.Z = Z; clear Z
        %PP_new.T = T; clear T
        %PP_new.acf = acf; clear acf
        %PP_new.bounds = bounds; clear bounds
        %PP_new.acf_over_thr = acf_over_thr; clear acf_over_thr
               
        % MONOVARIATE CENSORED (prof)
%         PP.Mu = Mu(:,1:UndSampFact:end); clear Mu 
%         PP.Thetap = Thetap(:,1:UndSampFact:end); clear Thetap
%         if opt.hasTheta0,PP.Theta0 = opt.Theta0; end
%         PP.meanRR = opt.meanRR(:,1:UndSampFact:end);
%         PP.Kappa = Kappa(:,1:UndSampFact:end); clear Kappa
%         PP.L = L(:,1:UndSampFact:end); clear L
%         PP.s2 = s2(:,1:UndSampFact:end); clear s2
%         PP.delta = delta*UndSampFact;
%         PP.t0 = opt.t0;
%         PP.LogLikel = opt.LogLikel(:,1:UndSampFact:end); clear opt
%         PP.KSdistance = KSdistance; clear KSdistance
%         PP.KSoutPerc = KSoutPerc; clear KSoutPercCov
%         PP.Z = Z; clear Z
%         PP.T = T; clear T
%         PP.acf = acf; clear acf
%         PP.bounds = bounds; clear bounds
%         PP.acf_over_thr = acf_over_thr; clear acf_over_thr
%         
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
% catch ME
%     disp('Failed')
%     disp(ME.message)
%     PP = deal([]);
%     PP.error = ME;
% end

end