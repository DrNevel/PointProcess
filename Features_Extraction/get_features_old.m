function [Features,RR_correction,Coherences] = get_features_old(selected_wf, varargin)
% [Features,RR_correction,Coherences] = get_features(selected_wf, 'fs',Fs,'get_PP',1,'get_TimeDomain',1)


% Last Updates:
% 25/03/2020 - Added Signals Resampling for Classical Features - MAX
% 15/05/2020 - Monovariate and Bivariate Models and Features Uploads - MAX
% 04/08/2020 - Addedd LF,HF Time-Varying Frequencies for Bivariate Spectra
% 24/12/2020 - Added cov_type='pat'

% Copyright (C) Maximiliano Mollura Riccardo Barbieri, 2019-2020.
% All Rights Reserved. See LICENSE.TXT for license details.
% maximiliano.mollura@polimi.it
% riccardo.barbieri@polimi.it



fs = 125;
toplot = 0;
get_PointProcess = 0;
get_PPfeat = 0;
% get_LSTM = 0;
get_TimeDomain = 0;
get_Stats = 0;
get_Freqs = 0;
get_Compl = 0;
get_mono = 0;
get_biva = 0;
get_multi = 0;
get_spectra = 0;
RR_correction =[];
Coherences = [];
PPorderMono = 9;
PPorderBiva = 3;
PPorderMulti = 3;
UndSampl = 1;
cov_type = 'syst';
sel_EKGR = 1;
rsmpl_ts = 0;
correct_diast=1;

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
        case 'get_Compl'
            get_Compl=varargin{i+1};
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
        case 'UndSampl'
            UndSampl = varargin{i+1};
        case 'cov_type'
            cov_type = varargin{i+1};
        case 'sel_EKGR'
            sel_EKGR = varargin{i+1};
        case 'rsmpl_ts'
            rsmpl_ts = varargin{i+1};
        case 'fres'
            fres = varargin{i+1};
        case 'correct_diast'
            correct_diast=varargin{i+1};
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
    end
    EKGR = S.annotations(sel_EKGR,:)./fs;
    ONSET = S.annotations(2,:)./fs;
    SYST = S.annotations(3,:)./fs;
    SYST_VAL = S.annotations(4,:);
    DIAST = S.annotations(5,:)./fs;
    DIAST_VAL = S.annotations(6,:);
    PTT = S.annotations(2,:)./fs - S.annotations(1,:)./fs;
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
catch
    error('Prepare Properly The Data!')
end

% Look at the drift in PTT (PPG)
% .....
%%%%%%%%%%%%%%%%%%%%%%%%%


EKGR_old = EKGR;
SYST_old = SYST;
DIAST_old = DIAST;
% Close holes in the Tachogram
[EKGR,SYST,DIAST,ONSET,RR_correction] = closeECG(EKGR,SYST,DIAST,ONSET);
TACO = diff(EKGR);

if rsmpl_ts
    TACO_ts = timeseries(TACO,EKGR(2:end));
    SYST_ts = timeseries(SYST_VAL,SYST);
    DIAST_ts = timeseries(DIAST_VAL,DIAST);
    PP_ts = timeseries(PP,DIAST);
    PTT_ts = timeseries(PTT,SYST);
    resfreq = 1/mean(TACO);
    % Resampling Using Re-sampled R-peaks Times
    taco_times = EKGR(2):resfreq:EKGR(end);
    syst_times = SYST(1):resfreq:SYST(end);
    diast_times =  DIAST(1):resfreq:DIAST(end);
    TACO_ts = resample(TACO_ts,taco_times,'linear');
    SYST_ts = resample(SYST_ts,syst_times,'linear');
    DIAST_ts = resample(DIAST_ts,diast_times,'linear');
    PTT_ts = resample(PTT_ts,syst_times,'linear');
    PP_ts = resample(PP_ts,diast_times,'linear');
    TACO_res = squeeze(TACO_ts.Data)';
    SYST_res = squeeze(SYST_ts.Data)';
    DIAST_res = squeeze(DIAST_ts.Data)';
    PP_res = squeeze(PP_ts.Data)';
    PTT_res = squeeze(PTT_ts.Data)';
else
    taco_times = EKGR(2:end);
    syst_times = SYST;
    diast_times =  DIAST;
    TACO_res = TACO;
    SYST_res = SYST_VAL;
    DIAST_res = DIAST_VAL;
    PP_res = PP;
    PTT_res = PTT;
end
% %% TO VISUALIZE
% plot(EKGR(2:end),TACO),hold on, plot(taco_times,TACO_res)
% plot(SYST,SYST_VAL),hold on, plot(syst_times,SYST_res)
% %%

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
end

% Extract Features
% PointProcess
if get_PointProcess
    if get_mono
        Features.PP.Mono = get_PP(EKGR, 'get_mono',1,'rr_corr', RR_correction,'get_spectra',get_spectra,'order', PPorderMono,'UndSampl',UndSampl);
        if get_spectra
            [Features.PP.Mono.RR.powLF, Features.PP.Mono.RR.powHF, Features.PP.Mono.RR.bal,...
                warn, Features.PP.Mono.RR.powVLF, Features.PP.Mono.RR.powTot,Features.PP.Mono.RR.S_RR,...
                Features.PP.Mono.RR.FS] =...
                hrv_indices(Features.PP.Mono.Thetap, Features.PP.Mono.s2, 1./Features.PP.Mono.meanRR, fres);
        end
    end
    if get_biva
        if get_mono
            Features.PP.Cov = get_PP(EKGR, 'get_biva',1,'cov_t',COV_T(:)','cov_v', COV_V(:)',...
                'rr_corr', RR_correction,'get_spectra',get_spectra,'order', PPorderBiva,'UndSampl',UndSampl,'Pmono',PPorderMono);
        else
            Features.PP.Cov = get_PP(EKGR, 'get_biva',1,'cov_t',COV_T(:)','cov_v', COV_V(:)',...
                'rr_corr', RR_correction,'get_spectra',get_spectra,'order', PPorderMulti,'UndSampl',UndSampl);
        end
        %         Features = get_PP(EKGR, 'get_biva',1,'cov_t',DIAST,'cov_v', DIAST_VAL,'rr_corr', RR_correction);
        %         Features = get_PP(EKGR, 'get_biva',1,'cov_t',SYST,'cov_v', PTT,'rr_corr', RR_correction);
        %         Features = get_PP(EKGR, 'get_biva',1,'cov_t',DIAST,'cov_v', PP,'rr_corr', RR_correction);
        
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
                end
                % Testing Plot
                %                 testplot(Features,EKGR_old, COV_V)
            end
        end
    end
    if get_multi
        Features.PP.Cov = get_PP(EKGR, 'get_multi',1,'cov_t',COV_T,'cov_v', COV_V,...
                'rr_corr', RR_correction,'get_spectra',get_spectra,'order', PPorderBiva,'UndSampl',UndSampl);
        %         Features = get_PP(EKGR, 'get_mono',0,'get_multi',1,'cov_t',{SYST,SYST},'cov_v', {SYST_VAL,PTT},'rr_corr', RR_correction);
    end
end

Features.sig_duration = DIAST(end)-EKGR(1); % seconds

% Time Domain
if get_TimeDomain
    [Features.CLASS.AVNN] = get_AVNN(TACO_res);
    [Features.CLASS.SDANN] = get_SDANN(TACO_res);
    [Features.CLASS.SDNN] = get_SDNN(TACO_res);
    [Features.CLASS.SDNNIDX] = get_SDNNIDX(TACO_res);
    [Features.CLASS.NN20, Features.CLASS.pNN20] = get_NN20(TACO_res);
    [Features.CLASS.NN50, Features.CLASS.pNN50] = get_NN50(TACO_res);
    [Features.CLASS.RMSSD] = get_RMSSD(TACO_res); % Output in mseconds
    [Features.CLASS.logRMSSD] = log(Features.CLASS.RMSSD);
    [Features.CLASS.SDSD] = get_SDSD(TACO_res);
    [Features.CLASS.SD1, Features.CLASS.SD2, Features.CLASS.SD_ratio, Features.CLASS.SD_prod]=get_poincare(TACO_res);
    [Features.CLASS.TRI,Features.CLASS.TINN] = get_TriangVals(TACO_res); % Rivedere  TINN
    [Features.CLASS.AVSS] = mean(SYST_res,'omitnan');
    [Features.CLASS.SDSS] = std(SYST_res,'omitnan');
    [Features.CLASS.AVDD] = mean(DIAST_res,'omitnan');
    [Features.CLASS.SDDD] = std(DIAST_res,'omitnan');
    [Features.CLASS.AVPP] = mean(PP_res,'omitnan');
    [Features.CLASS.SDPP] = std(PP_res,'omitnan');
    [Features.CLASS.AVPTT] = mean(PTT_res,'omitnan');
    [Features.CLASS.SDPTT] = std(PTT_res,'omitnan');
end

% Spectral Features
if get_Freqs
    % RR - HRV parameters
    [Features.CLASS.RR_TOTPWR,~,Features.CLASS.RR_VLF, Features.CLASS.RR_LF, Features.CLASS.RR_HF, ...
        Features.CLASS.RR_LFtoHF, Features.CLASS.RR_LFnu, Features.CLASS.RR_HFnu, Features.CLASS.RR_spect_slope,out] = get_spect2(TACO_res,taco_times,'pyulear');
    % SS - BPV parameters (Systolic Variability)
    [Features.CLASS.SS_TOTPWR,~,Features.CLASS.SS_VLF, Features.CLASS.SS_LF, Features.CLASS.SS_HF, ...
        ~, ~, ~, Features.CLASS.SS_spect_slope] = get_spect2(SYST_res,syst_times,'pyulear');
    % DD - BPV parameters (Diastolic Variability)
    [Features.CLASS.DD_TOTPWR,~,Features.CLASS.DD_VLF, Features.CLASS.DD_LF, Features.CLASS.DD_HF, ...
        ~, ~, ~, Features.CLASS.DD_spect_slope] = get_spect2(DIAST_res,diast_times,'pyulear');
    % PP  - BPV parameters (Pulse Pressure Variability)
    [Features.CLASS.PP_TOTPWR,~,Features.CLASS.PP_VLF, Features.CLASS.PP_LF, Features.CLASS.PP_HF, ...
        ~, ~, ~, Features.CLASS.PP_spect_slope] = get_spect2(PP_res,DIAST,'pyulear');
    % PTT - PTTVar parameters (PTT Variability)
    [Features.CLASS.PTT_TOTPWR,~,Features.CLASS.PTT_VLF, Features.CLASS.PTT_LF, Features.CLASS.PTT_HF, ...
        ~, ~, ~, Features.CLASS.PTT_spect_slope] = get_spect2(PTT_res,SYST,'pyulear');
end

% Complexity Measures
if get_Compl
    % Here use the non Interpolated TACOGRAM
    embed_dimension = 2;
    lag = 1;
    radius = 0.2*sqrt(trace(cov(TACO)));
    % RR
    if length(TACO) > 10000
        %         ratio = floor(length(RR)/10000);
        ratio = 10;
        [~,Features.CLASS.Alpha1,Features.CLASS.H,Features.CLASS.ApEn,Features.CLASS.sampen...
            ,Features.CLASS.CorrDim,Features.CLASS.LyapExp] = get_complexity(TACO(1:ratio:end), embed_dimension, lag, radius);
    else
        [~,Features.CLASS.Alpha1,Features.CLASS.H,Features.CLASS.ApEn,Features.CLASS.SampEn...
            ,Features.CLASS.CorrDim,Features.CLASS.LyapExp] = get_complexity(TACO, embed_dimension, lag, radius);
    end
    
end

Features.SIGS.ECG = EKGR;
Features.SIGS.SYST_VAL = SYST_VAL;
Features.SIGS.DIAST_VAL = DIAST_VAL;
Features.SIGS.SYST_T = SYST;
Features.SIGS.DIAST_T = DIAST;
Features.SIGS.PTT = PTT;
Features.SIGS.PP = PP;
Features.SIGS.MAP = MAP;
if rsmpl_ts
    Features.SIGS.Resampled.SYST_ts = SYST_ts;
    Features.SIGS.Resampled.DIAST_ts = DIAST_ts;
    Features.SIGS.Resampled.TACO_ts = TACO_ts;
    Features.SIGS.Resampled.PTT_ts = PTT_ts;
    Features.SIGS.Resampled.PP_ts = PP_ts;
end
Features.SIGS.ECG = EKGR;
Features.SIGS.SYST_VAL = SYST_VAL;
Features.SIGS.DIAST_VAL = DIAST_VAL;
Features.SIGS.SYST_T = SYST;
Features.SIGS.DIAST_T = DIAST;
Features.SIGS.PTT = PTT;
Features.SIGS.PP = PP;
Features.SIGS.MAP = MAP;
Features.SIGS.oldT.ECG = EKGR_old;
Features.SIGS.oldT.SYST_T = SYST_old;
Features.SIGS.oldT.DIAST_T = DIAST_old;





