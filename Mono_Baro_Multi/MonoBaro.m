function [outs,vout] = MonoBaro(RR_Theta, ord, s2_11, fsamp, prec, ret, lf_hf)
% [coher,fs,RSA_LF,RSA_HF,FRSA,PRSA_LF,PRSA_HF,RSA_all] = Baro(RR_Theta,COV_Theta, ord, s2_11, s2_22, fsamp, prec, modes)
% See R. Barbieri, G. Parati and J. P. Saul,
% "Closed- versus open-loop assessment of heart rate baroreflex,"
% in IEEE Engineering in Medicine and Biology Magazine, vol. 20, no. 2,
% pp. 33-42, March-April 2001.
%
%outs = {'GAIN_21_LF','GAIN_21_HF','PHASE_21_LF','PHASE_21_HF','GAIN_12_LF','GAIN_12_HF',...
%     'PHASE_12_LF','PHASE_12_HF','GCI_21_LF','GCI_21_HF','GCI_12_LF','GCI_12_HF',...
%     'COH_LF','COH_HF','F_LF','F_HF','S_RR_LF','S_RR_HF','S_BP_LF','S_RR_HF','S_CR_LF','S_CR_HF',...
%     'S_RR','S_BP','S_CR','COH','GAIN_2to1','PHASE_2to1',...
%     'GAIN_1to2','PHASE_1to2','FS','r_LF','r_HF','GCI_2to1','GCI_1to2','H11','H12','H21','H22','S_CR_conj'};
%
% Last Update: 15/5/2020
% Last Update: 04/8/2020 - Added SS_RR_LF,SS_RR_HF,SS_BP_LF,SS_BP_HF,SS_CR_LF,SS_RR_HF
% Last Update: 22/03/2021 - Normalize by the time-varying deltaF the spectra, *df

% Copyright (C) Maximiliano Mollura and Riccardo Barbieri, 2019-2020.
% All Rights Reserved. See LICENSE.TXT for license details.
% maximiliano.mollura@polimi.it
% riccardo.barbieri@polimi.it

% % conventional boundaries for LF and HF are assigned if not provided
% % customly
if ~exist('lf_hf','var')
    lf_hf = [.04 .15; .15 .45];
elseif isempty(lf_hf)
    lf_hf = [.04 .15; .15 .45];
end

% Checks
% % (code fails if parameters are not of the same order and same
% % "sampling")
[lt, N] = size(RR_Theta);
assert(lt == ord); % no n*ord as the model is monovariate, so n=1


% Freqeuncy Axis
FS = linspace(0, 1/2, prec); % normalized fsamp = from 0:1/2 %% prec ('fres' by backtracking, is the number of freq. partitions, 512 by default)
FS = FS'.*fsamp(:)'; % % N.W.: here fsamp is treated as a vector, but
                     % % backtracking in the code reveals that fsamp is
                     % % "1./Features.PP.Cov.meanRR" (scalar or vector?)
% % it would be nice to have a scalar 'fsamp', used as a scaling factor for
% % normalized frequencies, but the nature of "1./Features.PP.Cov.meanRR"
% % has to be understood

% % fsamp is later used in the function length(), supporting that
% % "1./Features.PP.Cov.meanRR" is indeed a vector

% % WHY IS "1./Features.PP.Cov.meanRR" A VECTOR? ISN'T fsamp SUPPOSED TO BE
% % SCALAR?

% variances - cross-variances assumed to be null

% a_0 set to be 0s
RR_Theta = [zeros(1,length(RR_Theta(1,:))); RR_Theta];

% ARs Coefficient
A = fft(RR_Theta(1:ord+1,:),prec*2,1);      % A11 is derived as the DFT of the TV RR AR parameters in the first eq. (considered as the first half of the RR param. vector)

% % to do the fft we need to do it up to f_samp (prec*2), but then we use
% % the transform up to f_nyquist (prec) 
A = A(1:prec,:);

% % Noise TFs are derived according to the theory
H = 1./(1 - A);

% Compute Spectra
% % noise (residual?) variances are derived as features from the PP model
S_RR = (abs(H).^2 .* s2_11)./fsamp;
% S_RR = ((conj(H).*H) .* s2_11)./fsamp;

% Compute Coherence
% COH = sqrt(S_CR.^2./(S_RR.*S_BP));

% Compute Baroreflex and FeedForward Gains
% % Gains are corrected to take into account "self-interference" (non sono
% sicuro del termine, ma "toglie" il contributo del segnale "self" al gain)
% GAIN_2to1 = sqrt( (S_RR - abs(H11).^2.*s2_11)./(S_BP-abs(H21).^2.*s2_11) ); % SBP to RR
% GAIN_1to2 = sqrt( (S_BP - abs(H22).^2.*s2_22)./(S_RR-abs(H12).^2.*s2_22) ); % RR to SBP
% PHASE_2to1 = angle(H12./H22);
% PHASE_1to2 = angle(H21./H11);
% 
% % Granger Causality Index
% % % monovariate parametric spectra used to derive GCI
% if exist('S_RR_mono','var')&&exist('S_BP_mono','var')
%     GCI_2to1 = log( (S_RR_mono./(S_RR-abs(H12).^2.*s2_22)) );
%     GCI_1to2 = log( (S_BP_mono./(S_BP-abs(H21).^2.*s2_11)) );
% elseif exist('S_RR_mono','var')&&~exist('S_BP_mono','var')
%     GCI_2to1 = log( (S_RR_mono./(S_RR-abs(H12).^2.*s2_22)) );
%     GCI_1to2 = log( (S_BP./(S_BP-abs(H21).^2.*s2_11)) );
% else
%     GCI_2to1 = log( (S_RR./(S_RR-abs(H12).^2.*s2_22)) );
%     GCI_1to2 = log( (S_BP./(S_BP-abs(H21).^2.*s2_11)) );
% end
% Frequency Grid
% UGUALIIIIIII!!!!
r_LF = lf_hf(1) < FS & FS < lf_hf(2); % % binary vector, true at LF indexes
r_HF = lf_hf(3) < FS & FS < lf_hf(4); % % binary vector, true at HF indexes
r_VLF = 0.003 < FS & FS <= lf_hf(1);   % % binary vector, true at VLF indexes

S_RR_TOT = NaN(1,length(fsamp));
S_RR_VLF = NaN(1,length(fsamp));
S_RR_LF = NaN(1,length(fsamp));
S_RR_HF = NaN(1,length(fsamp));

        for f = length(fsamp):-1:1
            df=(fsamp(f)/2)/prec;
            % Spectra as in mode 1
            if ~isnan(fsamp(f))
                S_RR_TOT(f) = sum(S_RR(:,f).*df,'omitnan');
                S_RR_VLF(f) = sum(S_RR(r_VLF(:,f),f).*df,'omitnan');
                S_RR_LF(f) = sum(S_RR(r_LF(:,f),f).*df,'omitnan');
                S_RR_HF(f) = sum(S_RR(r_HF(:,f),f).*df,'omitnan');
            end
        end
%     else
%         [GAIN_21_LF,GAIN_21_HF,GAIN_12_LF,GAIN_12_HF,PHASE_21_LF,PHASE_21_HF,PHASE_12_LF,...
%             PHASE_12_HF,GCI_21_LF,GCI_21_HF,GCI_12_LF,GCI_12_HF,COH_LF,COH_HF,F_LF,F_HF,S_RR_TOT,S_RR_VLF,S_RR_LF,...
%             S_RR_HF,S_BP_TOT,S_BP_VLF,S_BP_LF,S_BP_HF,S_CR_TOT,S_CR_VLF,S_CR_LF,S_CR_HF] = deal(0);
%     end
% else
%     [GAIN_21_LF,GAIN_21_HF,GAIN_12_LF,GAIN_12_HF,PHASE_21_LF,PHASE_21_HF,PHASE_12_LF,...
%             PHASE_12_HF,GCI_21_LF,GCI_21_HF,GCI_12_LF,GCI_12_HF,COH_LF,COH_HF,F_LF,F_HF,S_RR_TOT,S_RR_VLF,S_RR_LF,...
%             S_RR_HF,S_BP_TOT,S_BP_VLF,S_BP_LF,S_BP_HF,S_CR_TOT,S_CR_VLF,S_CR_LF,S_CR_HF] = deal(0);
%     
% end
outs = {'FS','r_LF','r_HF','H','S_RR','S_RR_TOT','S_RR_VLF','S_RR_LF','S_RR_HF'};

for out=1:length(ret)
    vout{out}=eval(outs{ret(out)});
end

% CONDENSED.GAIN_21_LF = GAIN_21_LF; clear GAIN_21_LF
% CONDENSED.GAIN_21_HF = GAIN_21_HF; clear GAIN_21_HF
% CONDENSED.GAIN_12_LF = GAIN_12_LF; clear GAIN_12_LF
% CONDENSED.GAIN_12_HF = GAIN_12_HF; clear GAIN_12_HF
% CONDENSED.GCI_21_LF = GCI_21_LF; clear GCI_21_LF
% CONDENSED.GCI_21_HF = GCI_21_HF; clear GCI_21_HF
% CONDENSED.GCI_12_LF = GCI_12_LF; clear GCI_12_LF
% CONDENSED.GCI_12_HF = GCI_12_HF; clear GCI_12_HF
% CONDENSED.F_LF = F_LF; clear F_LF
% CONDENSED.F_HF = F_HF; clear F_HF
%
% SPECTRA.S_RR = S_RR; clear S_RR
% SPECTRA.S_BP = S_BP; clear S_BP
% SPECTRA.S_CR = S_CR; clear S_CR
%
% COHERENCES.COH = COH; clear COH
%
% GAINS.GAIN_2to1 = GAIN_2to1; clear GAIN_2to1
% GAINS.GAIN_1to2 = GAIN_1to2; clear GAIN_1to2
%
% GCIS.GCI_2to1 = GCI_2to1; clear GCI_2to1
% GCIS.GCI_1to2 = GCI_1to2; clear GCI_1to2
%
% FREQ.FS = FS; clear FS
% FREQ.FS_LF = r_LF; clear r_LF
% FREQ.FS_HF = r_HF; clear r_HF
