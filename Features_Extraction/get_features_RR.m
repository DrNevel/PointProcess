function [Features,RR_correction,Coherences] = get_features_RR(selected_wf, varargin)
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
% 07/10/2021 - Only RR PP Estimation

% Copyright (C) Maximiliano Mollura Riccardo Barbieri, 2019-2020.
% All Rights Reserved. See LICENSE.TXT for license details.
% maximiliano.mollura@polimi.it
% riccardo.barbieri@polimi.it

%%%%

fs = 125;
toplot = 0;
get_PointProcess = 0;
get_TimeDomain = 0;
get_Freqs = 0;
mode='pyulear';
whole='';
get_Compl = 0;
get_spectra = 0;
RR_correction =[];
Coherences = [];
PPorderMono = 9;
W=90;
UndSampl = 1;
sel_EKGR = 1;
rsmpl_ts = 0;
fres=512;
every5min=1;

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
        case 'get_spectra'
            get_spectra=varargin{i+1};
        case 'PPorderMono'
            PPorderMono=varargin{i+1};
        case 'W'
            W=varargin{i+1};
        case 'UndSampl'
            UndSampl = varargin{i+1};
        case 'sel_EKGR'
            sel_EKGR = varargin{i+1};
            if sum(sel_EKGR == [3,5,7,9])
                warning('!!!!! Suspected sel_EKGR !!!!')
            end
        case 'rsmpl_ts'
            rsmpl_ts = varargin{i+1};
        case 'fres'
            fres = varargin{i+1};
        case 'every5min'
            every5min=varargin{i+1};
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
    elseif size(S.annotations,1)>=7
        EKGR = S.annotations(sel_EKGR,:)./fs; % from indexes S.annot to times EKGR
    end

catch
    error('Prepare Properly The Data!')
end

EKGR_old = EKGR;

% Close holes in the Tachogram
[EKGR,~,~,~,RR_correction] = closeECG(EKGR);
TACO = diff(EKGR);

if rsmpl_ts
    TACO_ts = timeseries(TACO,EKGR(2:end));
    resfreq = 1/mean(TACO);
    % Resampling Using Re-sampled R-peaks Times
    taco_times = EKGR(2):resfreq:EKGR(end);
    TACO_ts = resample(TACO_ts,taco_times,'linear');
    TACO_res = squeeze(TACO_ts.Data)';
else
    taco_times = EKGR(2:end);
    TACO_res = TACO;
end
% %% TO VISUALIZE
% plot(EKGR(2:end),TACO),hold on, plot(taco_times,TACO_res) 
% plot(SYST,SYST_VAL),hold on, plot(syst_times,SYST_res)
% %%

% Extract Features
% PointProcess
if get_PointProcess
        [Features.PP_new.Mono, Features.PP_U.Mono, Features.PP.Mono] = get_PP_edit(EKGR,'W',W,'get_mono',1,'rr_corr', RR_correction,'get_spectra',get_spectra,'order', PPorderMono,'UndSampl',UndSampl);
        
        if get_spectra
            % SPECTRA NEW (nostra)
            [Features.PP_new.Mono.RR.powLF, Features.PP_new.Mono.RR.powHF, Features.PP_new.Mono.RR.bal, ...
                warn, Features.PP_new.Mono.RR.powVLF, Features.PP_new.Mono.RR.powTot,...
                Features.PP_new.Mono.RR.S_RR, Features.PP_new.Mono.RR.FS] = ...
                hrv_indices(Features.PP_new.Mono.Thetap, Features.PP_new.Mono.s2, 1./Features.PP_new.Mono.meanRR, fres);
           
            % SPECTRA OLD UNCENSORED (prof modificata)
            [Features.PP_U.Mono.RR.powLF, Features.PP_U.Mono.RR.powHF, Features.PP_U.Mono.RR.bal, ...
                warn, Features.PP_U.Mono.RR.powVLF, Features.PP_U.Mono.RR.powTot,...
                Features.PP_U.Mono.RR.S_RR, Features.PP_U.Mono.RR.FS] = ...
                hrv_indices(Features.PP_U.Mono.Thetap, Features.PP_U.Mono.s2, 1./Features.PP_U.Mono.meanRR, fres);
            
            % SPECTRA CENSORED (prof)
            [Features.PP.Mono.RR.powLF, Features.PP.Mono.RR.powHF, Features.PP.Mono.RR.bal,...
                warn, Features.PP.Mono.RR.powVLF, Features.PP.Mono.RR.powTot,Features.PP.Mono.RR.S_RR,...
                Features.PP.Mono.RR.FS] =...
                hrv_indices(Features.PP.Mono.Thetap, Features.PP.Mono.s2, 1./Features.PP.Mono.meanRR, fres);
        end
        %         Features = get_PP(EKGR, 'get_biva',1,'cov_t',DIAST,'cov_v', DIAST_VAL,'rr_corr', RR_correction);
        %         Features = get_PP(EKGR, 'get_biva',1,'cov_t',SYST,'cov_v', PTT,'rr_corr', RR_correction);
        %         Features = get_PP(EKGR, 'get_biva',1,'cov_t',DIAST,'cov_v', PP,'rr_corr', RR_correction);
end

Features.sig_duration = EKGR(end)-EKGR(1); % seconds

% Classical Time Domain Features
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


Features.nbeats = length(EKGR);
Features.SIGS.ECG = EKGR;
if rsmpl_ts
    Features.SIGS.Resampled.TACO_ts = TACO_ts;
end
Features.SIGS.ECG = EKGR;
Features.SIGS.oldT.ECG = EKGR_old;


