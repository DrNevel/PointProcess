function [Features,RR_correction,Coherences] = get_features(selected_wf, varargin)
% [Features,RR_correction,Coherences] = get_features(selected_wf, 'fs',Fs,'get_PP',1,'get_TimeDomain',1)


% Last Updates:
% 25/03/2020 - Added Signals Resampling for Classical Features - MAX
% 15/05/2020 - Monovariate and Bivariate Models and Features Uploads - MAX
% 04/08/2020 - Addedd LF,HF Time-Varying Frequencies for Bivariate Spectra
% 24/12/2020 - Added cov_type='pat'
% 12/02/2021 - Changed PTT to PAT, corrected PAT changes after hole
%              Added SS,DD (as PRV) and SAP,DAP as BPV measures, changed
%              LF,HF, etc scales for HRV and PRV to msec^2
% 10/03/2021 - Updated LF,HF,... power integrals

% Copyright (C) Maximiliano Mollura Riccardo Barbieri, 2019-2020.
% All Rights Reserved. See LICENSE.TXT for license details.
% maximiliano.mollura@polimi.it
% riccardo.barbieri@polimi.it



fs = 125;
toplot = 0;
get_PointProcess = 0;
get_PPfeat = 0;
get_TimeDomain = 0;
get_Stats = 0;
get_Freqs = 0;
mode='pyulear';
whole='';
get_Compl = 0;
get_Compl_others = 0;
get_mono = 0;
get_biva = 0;
get_multi = 0;
get_spectra = 0;
RR_correction =[];
Coherences = [];
PPorderMono = 9;
PPorderBiva = 3;
PPorderMulti = 3;
W=90;
UndSampl = 1;
cov_type = 'syst';
sel_EKGR = 1;
rsmpl_ts = 0;
correct_diast=0;
fres=512;
compl_vals=[];
every5min=1;
hrv_pressure_variability=0;

for i = 1:2:length(varargin)
    switch varargin{i}
        case 'fs'
            fs=varargin{i+1};
        case 'toplot'
            toplot=varargin{i+1};
        case 'get_PP'
            get_PointProcess=varargin{i+1};
        case 'get_TimeDomain'
            get_TimeDomain=varargin{i+1};
        case 'get_Freqs'
            get_Freqs=varargin{i+1};
        case 'mode'
            mode=varargin{i+1};
        case 'whole'
            whole=varargin{i+1};
        case 'get_Compl'
            get_Compl=varargin{i+1};
        case 'get_Compl_others'
            get_Compl_others=varargin{i+1};
        case 'compl_vals'
            compl_vals=varargin{i+1};
        case 'get_mono'
            get_mono=varargin{i+1};
        case 'get_biva'
            get_biva=varargin{i+1};
        case 'get_multi'
            get_multi=varargin{i+1};
        case 'get_spectra'
            get_spectra=varargin{i+1};
        case 'PPorderMono'
            PPorderMono=varargin{i+1};
        case 'PPorderBiva'
            PPorderBiva=varargin{i+1};
        case 'PPorderMulti'
            PPorderMulti=varargin{i+1};
        case 'W'
            W=varargin{i+1};
        case 'UndSampl'
            UndSampl = varargin{i+1};
        case 'cov_type'
            cov_type = varargin{i+1};
        case 'sel_EKGR'
            sel_EKGR = varargin{i+1};
            if sum(sel_EKGR == [3,5,7,9])
                warning('!!!!! Suspected sel_EKGR !!!!')
            end
        case 'rsmpl_ts'
            rsmpl_ts = varargin{i+1};
        case 'fres'
            fres = varargin{i+1};
        case 'correct_diast'
            correct_diast=varargin{i+1};
        case 'every5min'
            every5min=varargin{i+1};
        case 'hrv_press'
            hrv_pressure_variability=varargin{i+1};
    end
end

%% START

% Data Extraction
try
    S = load(selected_wf);
    fn = fieldnames(S);
    if isstruct(S.(fn{1}))
        if isfield(S.(fn{1}),'annotations')
            S=S.(fn{1});
        end
    end
    if isfield(S,'Sig')
        S = S.Sig;
    end
    if isfield(S,'Fs')
        fs = S.Fs;
    elseif isfield(S,'fs')
        fs = S.fs;
    end
    if size(S.annotations,1)==1
        EKGR = S.annotations(1,:)./fs;        
    elseif size(S.annotations,1)==7 % Analysis with ECG-SAP
        EKGR = S.annotations(sel_EKGR,:)./fs;
        ONSET = S.annotations(2,:)./fs;
        ONSET_VAL = S.annotations(3,:);
        SYST = S.annotations(4,:)./fs;
        SYST_VAL = S.annotations(5,:);
        DIAST = S.annotations(6,:)./fs;
        DIAST_VAL = S.annotations(7,:);
        %PAT = S.annotations(2,:)./fs - S.annotations(1,:)./fs;
        PP = SYST_VAL-DIAST_VAL;
        MAP = (SYST_VAL+2*DIAST_VAL)./3;
%     elseif size(S.annotations,1)==?? % OTHER CONFIGURATIONS

    elseif size(S.annotations,1)>7 % Very Old Structure
        EKGR = S.annotations(sel_EKGR,:)./fs;
        ONSET = S.annotations(2,:)./fs;
        ONSET_VAL = [];%S.annotations(3,:);
        SYST = S.annotations(3,:)./fs;
        SYST_VAL = S.annotations(4,:);
        DIAST = S.annotations(5,:)./fs;
        DIAST_VAL = S.annotations(6,:);
        %PAT = S.annotations(2,:)./fs - S.annotations(1,:)./fs;
        PP = SYST_VAL-DIAST_VAL;

    if isfield(S,'abp')
        warning('MAP Momentaneously Changed!!!! Change for future studies')
        %         MAP = get_MAP60(S.abp,S.Fs);
        MAP = (SYST_VAL+2*DIAST_VAL)./3;
    else
        MAP = (SYST_VAL+2*DIAST_VAL)./3;
    end
    %         MAP = correct_abnormalities(MAP, 2, 3, 20*fs);%,10);
    %         All in seconds.
    end
catch
    error('Prepare Properly The Data!')
end

    

% Look at the drift in PTT (PPG)
% .....
%%%%%%%%%%%%%%%%%%%%%%%%%


EKGR_old = EKGR;
SYST_old = SYST;
DIAST_old = DIAST;
ONSET_old = ONSET;
% Close holes in the Tachogram
[EKGR,SYST,DIAST,ONSET,RR_correction] = closeECG(EKGR,SYST,DIAST,ONSET,correct_diast);
TACO = diff(EKGR);
PAT = ONSET-EKGR;

if rsmpl_ts
    TACO_ts = timeseries(TACO,EKGR(2:end));
    SYST_ts = timeseries(SYST_VAL,SYST);
    DIAST_ts = timeseries(DIAST_VAL,DIAST);
    PP_ts = timeseries(PP,DIAST);
    PAT_ts = timeseries(PAT,SYST);
    resfreq = 1/mean(TACO);
    % Resampling Using Re-sampled R-peaks Times
    taco_times = EKGR(2):resfreq:EKGR(end);
    onset_times = ONSET(1):resfreq:ONSET(end);
    syst_times = SYST(1):resfreq:SYST(end);
    diast_times =  DIAST(1):resfreq:DIAST(end);
    TACO_ts = resample(TACO_ts,taco_times,'linear');
    SYST_ts = resample(SYST_ts,syst_times,'linear');
    DIAST_ts = resample(DIAST_ts,diast_times,'linear');
    PAT_ts = resample(PAT_ts,syst_times,'linear');
    PP_ts = resample(PP_ts,diast_times,'linear');
    TACO_res = squeeze(TACO_ts.Data)';
    SYST_res = squeeze(SYST_ts.Data)';
    DIAST_res = squeeze(DIAST_ts.Data)';
    PP_res = squeeze(PP_ts.Data)';
    PAT_res = squeeze(PAT_ts.Data)';
    if ~isempty(ONSET_VAL)
        ONSET_ts = timeseries(ONSET_VAL,ONSET);
        ONSET_ts = resample(ONSET_ts,onset_times,'linear');
        ONSET_res = squeeze(ONSET_ts.Data)';
    else
        ONSET_ts = [];
        ONSET_res = [];
    end
else
    taco_times = EKGR(2:end);
    onset_times = ONSET;
    syst_times = SYST;
    diast_times =  DIAST;
    TACO_res = TACO;
    ONSET_res = ONSET_VAL;
    SYST_res = SYST_VAL;
    DIAST_res = DIAST_VAL;
    PP_res = PP;
    PAT_res = PAT;
end
% %% TO VISUALIZE
% plot(EKGR(2:end),TACO),hold on, plot(taco_times,TACO_res)
% plot(SYST,SYST_VAL),hold on, plot(syst_times,SYST_res)
% %%

% Select COVARIATE
if strcmp(cov_type,'syst')
    COV_T = SYST;
    COV_V = SYST_VAL;
elseif strcmp(cov_type,'diast')&&(correct_diast==1)
    EKGR = EKGR(2:end);
    COV_T = DIAST(1:end-1);
    COV_V = DIAST_VAL(1:end-1);
elseif strcmp(cov_type,'diast')&&(correct_diast==0)
    COV_T = DIAST;
    COV_V = DIAST_VAL;
elseif strcmp(cov_type,'pat')
    COV_T = ONSET;
    COV_V = ONSET-EKGR;
elseif strcmp(cov_type,'pat-syst')
    COV_T = [ONSET;SYST];
    COV_V = [ONSET-EKGR;SYST_VAL];
    get_mono=0;
    get_biva=0;
    get_multi=1;
elseif strcmp(cov_type,'rr_ss-ss')
    COV_T = [SYST(2:end);SYST(2:end)];
    COV_V = [diff(EKGR)-diff(SYST);SYST_VAL(2:end)]; %[RR-SS,SAP]
    EKGR = EKGR(2:end);
    get_mono=0;
    get_biva=0;
    get_multi=1;
end

%% Extract Features

% PointProcess
if get_PointProcess
    %%%%%%%%%%%%%%%% MONO PP %%%%%%%%%%%%%%%%
    if get_mono
        Features.PP.Mono = get_PP(EKGR,'W',W,'get_mono',1,'rr_corr', RR_correction,'get_spectra',get_spectra,'order', PPorderMono,'UndSampl',UndSampl);
        if get_spectra
            [Features.PP.Mono.RR.powLF, Features.PP.Mono.RR.powHF, Features.PP.Mono.RR.bal,...
                warn, Features.PP.Mono.RR.powVLF, Features.PP.Mono.RR.powTot,Features.PP.Mono.RR.S_RR,...
                Features.PP.Mono.RR.FS] =...
                hrv_indices(Features.PP.Mono.Thetap, Features.PP.Mono.s2, 1./Features.PP.Mono.meanRR, fres);
        end
    end
    %%%%%%%%%%%%%%%% BIVA PP %%%%%%%%%%%%%%%%
    if get_biva
        if get_mono 
            % (both mono and biva)
            Features.PP.Cov = get_PP(EKGR, 'W',W,'get_biva',1,'cov_t',COV_T(:)','cov_v', COV_V(:)',...
                'rr_corr', RR_correction,'get_spectra',get_spectra,'order', PPorderBiva,'UndSampl',UndSampl,'Pmono',PPorderMono);
        else        
            % (only biva)
            Features.PP.Cov = get_PP(EKGR, 'W',W,'get_biva',1,'cov_t',COV_T(:)','cov_v', COV_V(:)',...
                'rr_corr', RR_correction,'get_spectra',get_spectra,'order', PPorderBiva,'UndSampl',UndSampl);
        end        
        % Features = get_PP(EKGR, 'get_biva',1,'cov_t',DIAST,'cov_v', DIAST_VAL,'rr_corr', RR_correction);
        % Features = get_PP(EKGR, 'get_biva',1,'cov_t',SYST,'cov_v', PTT,'rr_corr', RR_correction);
        % Features = get_PP(EKGR, 'get_biva',1,'cov_t',DIAST,'cov_v', PP,'rr_corr', RR_correction);        
        if get_spectra
            if strcmp(cov_type,'syst')
                % Baroreflex - PP
                if get_mono
                    [Features.PP.Cov.SAP.powLF, Features.PP.Cov.SAP.powHF, Features.PP.Cov.SAP.bal,...
                        warn, Features.PP.Cov.SAP.powVLF, Features.PP.Cov.SAP.powTot,Features.PP.Cov.SAP.S_BP,...
                        Features.PP.Cov.SAP.FS] =...
                        hrv_indices(Features.PP.Cov.Thetap_cov, Features.PP.Cov.Var_cov, 1./Features.PP.Cov.meanRR, fres);
                    
                    [lbl,Features.PP.Cov.Baro.GAIN_21_LF,Features.PP.Cov.Baro.GAIN_21_HF,...
                        Features.PP.Cov.Baro.PHASE_21_LF,Features.PP.Cov.Baro.PHASE_21_HF,...
                        Features.PP.Cov.Baro.GAIN_12_LF,Features.PP.Cov.Baro.GAIN_12_HF,...
                        Features.PP.Cov.Baro.PHASE_12_LF,Features.PP.Cov.Baro.PHASE_12_HF,...
                        Features.PP.Cov.Baro.GCI_21_LF,Features.PP.Cov.Baro.GCI_21_HF,...
                        Features.PP.Cov.Baro.GCI_12_LF,Features.PP.Cov.Baro.GCI_12_HF,...
                        Features.PP.Cov.Baro.COH_LF,Features.PP.Cov.Baro.COH_HF,Features.PP.Cov.Baro.F_LF,...
                        Features.PP.Cov.Baro.F_HF,Features.PP.Cov.Baro.S_RR_TOT,Features.PP.Cov.Baro.S_RR_VLF,...
                        Features.PP.Cov.Baro.S_RR_LF,Features.PP.Cov.Baro.S_RR_HF,...
                        Features.PP.Cov.Baro.S_BP_TOT,Features.PP.Cov.Baro.S_BP_VLF,...
                        Features.PP.Cov.Baro.S_BP_LF,Features.PP.Cov.Baro.S_BP_HF,...
                        Features.PP.Cov.Baro.S_CR_TOT,Features.PP.Cov.Baro.S_CR_VLF,Features.PP.Cov.Baro.S_CR_LF,...
                        Features.PP.Cov.Baro.S_CR_HF,Features.PP.Cov.Baro.S_RR,Features.PP.Cov.Baro.S_BP,...
                        Features.PP.Cov.Baro.S_CR,Features.PP.Cov.Baro.COH,Features.PP.Cov.Baro.GAIN_2to1,...
                        Features.PP.Cov.Baro.PHASE_2to1,Features.PP.Cov.Baro.GAIN_1to2,Features.PP.Cov.Baro.PHASE_1to2,...
                        Features.PP.Cov.Baro.FS,Features.PP.Cov.Baro.r_LF,Features.PP.Cov.Baro.r_HF,...
                        Features.PP.Cov.Baro.GCI_2to1,Features.PP.Cov.Baro.GCI_1to2] = ...,H11,H12,H21,H22,S_CR_conj] = ...
                        ...
                        Baro_vect(Features.PP.Cov.Thetap,Features.PP.Cov.Thetap2(:,:,2),PPorderBiva,...
                        Features.PP.Cov.s2,Features.PP.Cov.cov_s2(2,:),1./Features.PP.Cov.meanRR,...
                        fres,1,1:46,Features.PP.Mono.RR.S_RR,Features.PP.Cov.SAP.S_BP);
                else
                    [lbl,Features.PP.Cov.Baro.GAIN_21_LF,Features.PP.Cov.Baro.GAIN_21_HF,...
                        Features.PP.Cov.Baro.PHASE_21_LF,Features.PP.Cov.Baro.PHASE_21_HF,...
                        Features.PP.Cov.Baro.GAIN_12_LF,Features.PP.Cov.Baro.GAIN_12_HF,...
                        Features.PP.Cov.Baro.PHASE_12_LF,Features.PP.Cov.Baro.PHASE_12_HF,...
                        Features.PP.Cov.Baro.GCI_21_LF,Features.PP.Cov.Baro.GCI_21_HF,...
                        Features.PP.Cov.Baro.GCI_12_LF,Features.PP.Cov.Baro.GCI_12_HF,...
                        Features.PP.Cov.Baro.COH_LF,Features.PP.Cov.Baro.COH_HF,Features.PP.Cov.Baro.F_LF,...
                        Features.PP.Cov.Baro.F_HF,Features.PP.Cov.Baro.S_RR_TOT,Features.PP.Cov.Baro.S_RR_VLF,...
                        Features.PP.Cov.Baro.S_RR_LF,Features.PP.Cov.Baro.S_RR_HF,...
                        Features.PP.Cov.Baro.S_BP_TOT,Features.PP.Cov.Baro.S_BP_VLF,...
                        Features.PP.Cov.Baro.S_BP_LF,Features.PP.Cov.Baro.S_BP_HF,...
                        Features.PP.Cov.Baro.S_CR_TOT,Features.PP.Cov.Baro.S_CR_VLF,Features.PP.Cov.Baro.S_CR_LF,...
                        Features.PP.Cov.Baro.S_CR_HF,Features.PP.Cov.Baro.S_RR,Features.PP.Cov.Baro.S_BP,...
                        Features.PP.Cov.Baro.S_CR,Features.PP.Cov.Baro.COH,Features.PP.Cov.Baro.GAIN_2to1,...
                        Features.PP.Cov.Baro.PHASE_2to1,Features.PP.Cov.Baro.GAIN_1to2,Features.PP.Cov.Baro.PHASE_1to2,...
                        Features.PP.Cov.Baro.FS,Features.PP.Cov.Baro.r_LF,Features.PP.Cov.Baro.r_HF,...
                        Features.PP.Cov.Baro.GCI_2to1,Features.PP.Cov.Baro.GCI_1to2] = ...,H11,H12,H21,H22,S_CR_conj] = ...
                        ...
                        Baro_vect(Features.PP.Cov.Thetap,Features.PP.Cov.Thetap2(:,:,2),PPorderBiva,...
                        Features.PP.Cov.s2,Features.PP.Cov.cov_s2(2,:),1./Features.PP.Cov.meanRR,...
                        fres,1,1:46);
                        %RR_Theta = Features.PP.Cov.Thetap
                        %COV_Theta = Features.PP.Cov.Thetap2(:,:,2)
                        %ord = PPorderBiva
                        %fsamp = 1./Features.PP.Cov.meanRR
                        %prec = fres
                        %ret = 1:46
                        %mode = 1
                end                
                % Testing Plot
                % testplot(Features,EKGR_old, COV_V)
            end
        end
    end
    
    %%%%%%%%%%%%%%%% MULTI PP %%%%%%%%%%%%%%%%
    if get_multi
        Features.PP.Cov = get_PP(EKGR, 'W',W,'get_multi',1,'cov_t',COV_T,'cov_v', COV_V,...
            'rr_corr', RR_correction,'get_spectra',get_spectra,'order', PPorderMulti,'UndSampl',UndSampl);
        % Features = get_PP(EKGR, 'get_mono',0,'get_multi',1,'cov_t',{SYST,SYST},'cov_v', {SYST_VAL,PTT},'rr_corr', RR_correction);
    
        %%%%%%%%%%%%%%%%% SPETTRI MULTIVAL %%%%%%%%%%%%%%%%%
        
    end
    
end

Features.sig_duration = DIAST(end)-EKGR(1); % seconds

% Time Domain
if get_TimeDomain
    [Features.CLASS.AVNN] = get_AVNN(TACO_res);
    [Features.CLASS.SDANN] = get_SDANN(TACO_res);
    [Features.CLASS.SDNN] = get_SDNN(TACO_res);
    [Features.CLASS.SDNNIDX] = get_SDNNIDX(TACO_res);
    [Features.CLASS.RMSSD] = get_RMSSD(TACO_res); % Output in mseconds
    [Features.CLASS.logRMSSD] = log(Features.CLASS.RMSSD);
    [Features.CLASS.SDSD] = get_SDSD(TACO_res);
    [Features.CLASS.SD1, Features.CLASS.SD2, Features.CLASS.SD_ratio, Features.CLASS.SD_prod]=get_poincare(TACO_res);
    [Features.CLASS.TRI,Features.CLASS.TINN] = get_TriangVals(TACO_res); % Rivedere  TINN
    [Features.CLASS.AVSAP] = mean(SYST_res,'omitnan');
    [Features.CLASS.SDSAP] = std(SYST_res,'omitnan');
    [Features.CLASS.AVDAP] = mean(DIAST_res,'omitnan');
    [Features.CLASS.SDDAP] = std(DIAST_res,'omitnan');
    [Features.CLASS.AVPP] = mean(PP_res,'omitnan');
    [Features.CLASS.SDPP] = std(PP_res,'omitnan');
    [Features.CLASS.AVPAT] = mean(PAT_res,'omitnan');
    [Features.CLASS.SDPAT] = std(PAT_res,'omitnan');
    
    Features.CLASS.NN20=[];
    Features.CLASS.NN50=[];
    Features.CLASS.pNN20=[];
    Features.CLASS.pNN50=[];
    if every5min
        cumul=0;
        start=1;
        for i=1:length(TACO_res)
            cumul=cumul+TACO_res(i);
            if (cumul>=300)||(cumul>=240&&i==length(TACO_res)) %admit last segment of 4 minutes
                ends=i;
                [NN20, pNN20] = get_NN20(TACO_res(start:ends));
                Features.CLASS.NN20=[Features.CLASS.NN20,NN20];
                Features.CLASS.pNN20=[Features.CLASS.pNN20,pNN20];
                [NN50, pNN50] = get_NN50(TACO_res(start:ends));
                Features.CLASS.NN50=[Features.CLASS.NN50,NN50];
                Features.CLASS.pNN50=[Features.CLASS.pNN50,pNN50];
                cumul=0;
                start=i+1;
            end
        end
        Features.CLASS.NN20=mean(Features.CLASS.NN20,'omitnan');
        Features.CLASS.pNN20=mean(Features.CLASS.pNN20,'omitnan');
        Features.CLASS.NN50=mean(Features.CLASS.NN50,'omitnan');
        Features.CLASS.pNN50=mean(Features.CLASS.pNN50,'omitnan');
    else
        [Features.CLASS.NN20, Features.CLASS.pNN20] = get_NN20(TACO_res);
        [Features.CLASS.NN50, Features.CLASS.pNN50] = get_NN50(TACO_res);
    end
end

% Spectral Features
if get_Freqs
    % RR - HRV parameters
    [Features.CLASS.RR_TOTPWR,~,Features.CLASS.RR_VLF, Features.CLASS.RR_LF, Features.CLASS.RR_HF, ...
        Features.CLASS.RR_LFtoHF, Features.CLASS.RR_LFnu, Features.CLASS.RR_HFnu, Features.CLASS.RR_spect_slope,Features.SPEC.RR] = get_spect2(TACO_res.*1000,taco_times,mode,whole);
    if hrv_pressure_variability
        % SS - RR parameters (Pulse Rate Variability)
        [Features.CLASS.SS_TOTPWR,~,Features.CLASS.SS_VLF, Features.CLASS.SS_LF, Features.CLASS.SS_HF, ...
            Features.CLASS.SS_LFtoHF, Features.CLASS.SS_LFnu, Features.CLASS.SS_HFnu, Features.CLASS.SS_spect_slope,Features.SPEC.SS] = get_spect2(diff(syst_times).*1000,syst_times(2:end),mode,whole);
        % OO - RR parameters (Pulse Rate Variability)
        [Features.CLASS.OO_TOTPWR,~,Features.CLASS.OO_VLF, Features.CLASS.OO_LF, Features.CLASS.OO_HF, ...
            Features.CLASS.OO_LFtoHF, Features.CLASS.OO_LFnu, Features.CLASS.OO_HFnu, Features.CLASS.OO_spect_slope,Features.SPEC.OO] = get_spect2(diff(onset_times).*1000,syst_times(2:end),mode,whole);
        % DD - RR parameters (Pulse Rate Variability)
        [Features.CLASS.DD_TOTPWR,~,Features.CLASS.DD_VLF, Features.CLASS.DD_LF, Features.CLASS.DD_HF, ...
            Features.CLASS.DD_LFtoHF, Features.CLASS.DD_LFnu, Features.CLASS.DD_HFnu, Features.CLASS.DD_spect_slope,Features.SPEC.DD] = get_spect2(diff(diast_times).*1000,diast_times(2:end),mode,whole);
    end
    % PAT - PATVar parameters (PAT Variability)
    [Features.CLASS.PAT_TOTPWR,~,Features.CLASS.PAT_VLF, Features.CLASS.PAT_LF, Features.CLASS.PAT_HF, ...
        ~, ~, ~, Features.CLASS.PAT_spect_slope,Features.SPEC.PAT] = get_spect2(PAT_res.*1000,SYST,mode,whole);
    % SAP - BPV parameters (Systolic Variability)
    [Features.CLASS.SAP_TOTPWR,~,Features.CLASS.SAP_VLF, Features.CLASS.SAP_LF, Features.CLASS.SAP_HF, ...
        ~, ~, ~, Features.CLASS.SAP_spect_slope,Features.SPEC.SAP] = get_spect2(SYST_res,syst_times,mode,whole);
    % ONP - RR parameters (Pulse Rate Variability)
    [Features.CLASS.ONP_TOTPWR,~,Features.CLASS.ONP_VLF, Features.CLASS.ONP_LF, Features.CLASS.ONP_HF, ...
        ~, ~, ~, Features.CLASS.ONP_spect_slope,Features.SPEC.ONP] = get_spect2(ONSET_res,onset_times,mode,whole);
    % DAP - BPV parameters (Diastolic Variability)
    [Features.CLASS.DAP_TOTPWR,~,Features.CLASS.DAP_VLF, Features.CLASS.DAP_LF, Features.CLASS.DAP_HF, ...
        ~, ~, ~, Features.CLASS.DAP_spect_slope,Features.SPEC.DAP] = get_spect2(DIAST_res,diast_times,mode,whole);
    % PP  - BPV parameters (Pulse Pressure Variability)
    [Features.CLASS.PP_TOTPWR,~,Features.CLASS.PP_VLF, Features.CLASS.PP_LF, Features.CLASS.PP_HF, ...
        ~, ~, ~, Features.CLASS.PP_spect_slope,Features.SPEC.PP] = get_spect2(PP_res,DIAST,mode,whole);
end

% Complexity Measures
if get_Compl
    % Here use the non Interpolated TACOGRAM
    embed_dimension = 2;
    lag = 1;
    radius = 0.2*sqrt(trace(cov(TACO)));
    % RR
    if length(TACO) > 10000
        ratio = 10;
    else
        ratio=1;
    end
    [~,Features.CLASS.Alpha1,Features.CLASS.H,Features.CLASS.ApEn,Features.CLASS.SampEn...
        ,Features.CLASS.CorrDim,Features.CLASS.LyapExp] = get_complexity(TACO(1:ratio:end), embed_dimension, lag, radius);
end

if get_Compl_others
    if length(TACO) > 10000
        ratio = 10;
    else
        ratio=1;
    end
    % Cross-Sample Entropy RR-SAP
    [Features.CLASS.XEntropy_RR_SAP,Features.XENTR.XEnts] = cross_sampen(TACO(1:ratio:end),SYST_VAL(2:ratio:end),3);
    [Features.CLASS.XEntropy_RR_SS,Features.XENTR.XEnts] = cross_sampen(TACO(1:ratio:end),diff(SYST(1:ratio:end)),3);
    [Features.CLASS.XEntropy_SSmRR,Features.XENTR.XEnts] = cross_sampen(TACO(1:ratio:end),diff(SYST(1:ratio:end))-TACO(1:ratio:end),3);
    % Sample Entropy on Other Signals
    if ~isempty(compl_vals)
        embed_dimensions=compl_vals.emb_dims_vect;
        lags=compl_vals.lags;
        radius = compl_vals.rads;
    else
        embed_dimensions=[2,3];
        lags=[1];
        radius=[0.2*sqrt(trace(cov(SYST_VAL)));...
            0.2*sqrt(trace(cov(DIAST_VAL)));...
            0.2*sqrt(trace(cov(PAT)))];
    end
    
    for embed_dimension=embed_dimensions
        row=embed_dimension-(min(embed_dimensions)-1);
        for lag=lags
            
            [~,Features.CLASS.SYST.Alpha1(row,lag),Features.CLASS.SYST.H(row,lag)...
                ,Features.CLASS.SYST.ApEn(row,lag),Features.CLASS.SYST.SampEn(row,lag)...
                ,Features.CLASS.SYST.CorrDim(row,lag),Features.CLASS.SYST.LyapExp(row,lag)...
                ] = get_complexity(SYST_VAL(1:ratio:end), embed_dimension, lag, radius(1));
            
            [~,Features.CLASS.DIAST.Alpha1(row,lag),Features.CLASS.DIAST.H(row,lag)...
                ,Features.CLASS.DIAST.ApEn(row,lag),Features.CLASS.DIAST.SampEn(row,lag)...
                ,Features.CLASS.DIAST.CorrDim(row,lag),Features.CLASS.DIAST.LyapExp(row,lag)...
                ] = get_complexity(DIAST_VAL(1:ratio:end), embed_dimension, lag, radius(2));
            
            [~,Features.CLASS.PAT.Alpha1(row,lag),Features.CLASS.PAT.H(row,lag)...
                ,Features.CLASS.PAT.ApEn(row,lag),Features.CLASS.PAT.SampEn(row,lag)...
                ,Features.CLASS.PAT.CorrDim(row,lag),Features.CLASS.PAT.LyapExp(row,lag)...
                ] = get_complexity(PAT(1:ratio:end), embed_dimension, lag, radius(3));
        end
    end
end

Features.nbeats = length(EKGR);
Features.SIGS.ECG = EKGR;
Features.SIGS.SYST_VAL = SYST_VAL;
Features.SIGS.DIAST_VAL = DIAST_VAL;
Features.SIGS.SYST_T = SYST;
Features.SIGS.DIAST_T = DIAST;
Features.SIGS.PAT = PAT;
Features.SIGS.PP = PP;
Features.SIGS.MAP = MAP;
if rsmpl_ts
    Features.SIGS.Resampled.SYST_ts = SYST_ts;
    Features.SIGS.Resampled.DIAST_ts = DIAST_ts;
    Features.SIGS.Resampled.TACO_ts = TACO_ts;
    Features.SIGS.Resampled.PAT_ts = PAT_ts;
    Features.SIGS.Resampled.PP_ts = PP_ts;
end
Features.SIGS.ECG = EKGR;
Features.SIGS.ONSET = ONSET;
Features.SIGS.ONSET_VAL = ONSET_VAL;
Features.SIGS.SYST_VAL = SYST_VAL;
Features.SIGS.DIAST_VAL = DIAST_VAL;
Features.SIGS.SYST_T = SYST;
Features.SIGS.DIAST_T = DIAST;
Features.SIGS.PAT = PAT;
Features.SIGS.PP = PP;
Features.SIGS.MAP = MAP;
Features.SIGS.oldT.ECG = EKGR_old;
Features.SIGS.oldT.ONSET_T = ONSET_old;
Features.SIGS.oldT.SYST_T = SYST_old;
Features.SIGS.oldT.DIAST_T = DIAST_old;





