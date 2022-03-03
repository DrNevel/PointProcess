function [outs,vout] = TriBaroShort(RR_Theta,COV_Theta, ord, s2_11, s2_22, s2_33, fsamp, prec, modes,ret,lf_hf,S_RR_mono,S_BP_mono,S_PAT_mono)
% [coher,fs,RSA_LF,RSA_HF,FRSA,PRSA_LF,PRSA_HF,RSA_all] = Baro(RR_Theta,COV_Theta, ord, s2_11, s2_22, fsamp, prec, modes)
% See R. Barbieri, G. Parati and J. P. Saul,
% "Closed- versus open-loop assessment of heart rate baroreflex,"
% in IEEE Engineering in Medicine and Biology Magazine, vol. 20, no. 2,
% pp. 33-42, March-April 2001.
%
% Let's begin with the indexes with clear and reported multivariate
% (trivariate) formulas, such as 'H' T.F.s, spectra, coherences:
% See M. Orini, G. Valenza, L. Citi and R. Barbieri,
% "Tetravariate Point-Process Model for the Continuous Characterization of
% Cardiovascular-Respiratory Dynamics during Passive Postural Changes"
% in Computing in Cardiology, 2012, pp. 273-276.
% 
%outs = {'GAIN_21_LF','GAIN_21_HF','PHASE_21_LF','PHASE_21_HF','GAIN_12_LF','GAIN_12_HF',...
%     'PHASE_12_LF','PHASE_12_HF','COH_LF','COH_HF','F_LF','F_HF',...
%     'S_RR_LF','S_RR_HF','S_BP_LF','S_RR_HF','S_CR_LF','S_CR_HF',...
%     'S_RR','S_BP','S_CR','COH','GAIN_2to1','PHASE_2to1',...
%     'GAIN_1to2','PHASE_1to2','FS','r_LF','r_HF','H11','H12','H21','H22','S_CR_conj'};
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
% TO CHECK: if the interpretation is correct, this comes from the
% monovariate PP optimization; the code has to be run in debug mode to
% check this fact, as it would be really helpful also in the SYMSIG part.
[lc, Nc, Ncov] = size(COV_Theta);
assert(Ncov == 3); % this code is thought explicitly to handle the trivariate problem
% TO CHECK: when data is available, check whether 'Features.PP.Cov.Thetap'
% comes from the monovariate PP optimization and 'Features.PP.Cov.Thetap2'
% comes from the bivariate one (in 'get_features.m'); from a quick look in
% 'get_PP.m' it seems to be the case, but some runs are needed.
assert(N == Nc);
assert(lt == lc);
% NW: probably 'Features.PP.Cov.Thetap' and 'Features.PP.Cov.Thetap2' are
% are respectively obtained from PP inv_gauss optim. and WLS gauss. estim.
assert(lt == 3*ord);


% Freqeuncy Axis
FS = linspace(0, 1/2, prec); % normalized fsamp = from 0:1/2 %% prec ('fres' by backtracking, is the number of freq. partitions, 512 by default)
FS = FS'.*fsamp(:)';

% % 'fsamp' is a vector of the time variant sampling frequency (RR-beats)
% % (1./Features.PP.Cov.meanRR); doubts that brought to this conclusion can
% % be seen in the comments in 'Baro.m'

% variances - cross-variances assumed to be null

% ARs Coefficient
COV_Theta1=COV_Theta(:,:,1);
COV_Theta2=COV_Theta(:,:,2);
COV_Theta3=COV_Theta(:,:,3);

index=find(isnan(RR_Theta(1,:))==0,1,'first')-1;

% a_0(s) to be changed into zeroes (to not get artifacts/shifts in spectra)
A0=zeros(1,length(RR_Theta(1,:)));
A0(1:index)=NaN;
RR_Theta = [A0;RR_Theta];
COV_Theta1 = [A0;COV_Theta(:,:,1)];
COV_Theta2 = [A0;COV_Theta(:,:,2)];
COV_Theta3 = [A0;COV_Theta(:,:,3)];

COV_Theta=cat(3, COV_Theta1, COV_Theta2, COV_Theta3);

% Aij are derived as the DFT of the TV RR AR parameters
% applied[if this does not work, use 'ord+1' and 'ord+2' as in the comment below]
A11 = fft(RR_Theta(1:ord+1,:),prec*2,1); % order shifted by 1 as a_0 is included
A12 = fft(RR_Theta(ord+2:2*ord+1,:),prec*2,1); % all indexes adequately shifted
A13 = fft(RR_Theta((2*ord+2):end,:),prec*2,1);
A21 = fft(COV_Theta(2:ord+1,:,2),prec*2,1);
A22 = fft([COV_Theta(1,:,2);COV_Theta(ord+2:2*ord+1,:,2)],prec*2,1);
A23 = fft(COV_Theta((2*ord+2):end,:,2),prec*2,1);
A31 = fft(COV_Theta(2:ord+1,:,3),prec*2,1);
A32 = fft(COV_Theta(ord+2:2*ord+1,:,3),prec*2,1);
A33 = fft([COV_Theta(1,:,3);COV_Theta((2*ord+2):end,:,3)],prec*2,1);
% [N.W.]: every Aij is zero padded up to 'prec*2', so the results will have
% the same dimensions even though Aii elements also consider a_0
% coefficients.

% A11 = fft(RR_Theta(1:P+1),prec*2,1);
% A12 = fft(RR_Theta(P+2:(2*P+1)),prec*2,1);
% A13 = fft(RR_Theta((2*P+2):end),prec*2,1);
% A21 = fft(COV_Theta(2:P+1),prec*2,1);
% A22 = fft([COV_Theta(1);COV_Theta(P+2:end)],prec*2,1);
% A23 = fft(COV_Theta((2*P+2):end),prec*2,1);
% A31 = fft(COV_Theta(2:P+1),prec*2,1);
% A32 = fft(COV_Theta(P+2:(2*P+1)),prec*2,1);
% A33 = fft([COV_Theta(1);COV_Theta((2*P+2):end)],prec*2,1);

% % to do the fft we need to do it up to f_samp (prec*2), but then we use
% % the transform up to f_nyquist (prec) 
A11 = A11(1:prec,:);
A12 = A12(1:prec,:);
A13 = A13(1:prec,:);
A21 = A21(1:prec,:);
A22 = A22(1:prec,:);
A23 = A23(1:prec,:);
A31 = A31(1:prec,:);
A32 = A32(1:prec,:);
A33 = A33(1:prec,:);

% % Noise TFs are derived according to the theory (bivariate)

% % Deductively, the trivariate forms should be the following:
% using subs(), the denominator of the following formulas and the
% determinant of A (refer to google doc 'tribaro') give the same result

% H is derived as inv(I - A) (in the case of errors, check TriBaro)
detA = (1-A11).*(1-A22).*(1-A33)-(A12.*A23.*A31)-(A13.*A21.*A32)-((1-A11).*A23.*A32)-((1-A33).*A12.*A21)-((1-A22).*A13.*A31);
H11  = ((A23.*A32)-((1-A22).*(1-A33)))./detA;
H12  = -((A12.*(1-A33))+(A13.*A32))./detA;
H13  = -((A13.*(1-A22))+(A12.*A23))./detA;
H21  = -((A21.*(1-A33))+(A23.*A31))./detA;
H22  = ((A13.*A31)-((1-A11).*(1-A33)))./detA;
H23  = -((A23.*(1-A11))+(A13.*A21))./detA;
H31  = -((A31.*(1-A22))+(A21.*A32))./detA;
H32  = -((A32.*(1-A11))+(A12.*A31))./detA;
H33  = ((A12.*A21)-((1-A11).*(1-A22)))./detA;

% Compute Spectra
% % noise (residual?) variances are derived as features from the PP model
% % Here should appear also a s2_33
S_RR  = abs(H11).^2 .* s2_11./fsamp + abs(H12).^2 .* s2_22./fsamp + abs(H13).^2 .* s2_33./fsamp;
S_BP  = abs(H21).^2 .* s2_11./fsamp + abs(H22).^2 .* s2_22./fsamp + abs(H23).^2 .* s2_33./fsamp;
S_PAT = abs(H31).^2 .* s2_11./fsamp + abs(H32).^2 .* s2_22./fsamp + abs(H33).^2 .* s2_33./fsamp;
% Controllare 'TriBaro_Redux.m' per vedere i primi risultati (erronei?)
% dei cross-spettri
S_CR12 = H21.*s2_11.*conj(H11) + H22.*s2_22.*conj(H21) + H23.*s2_33.*conj(H31);
S_CR13 = H31.*s2_11.*conj(H11) + H32.*s2_22.*conj(H21) + H33.*s2_33.*conj(H31);
S_CR23 = H31.*s2_11.*conj(H12) + H32.*s2_22.*conj(H22) + H33.*s2_33.*conj(H32);
S_CR12_conj = H11.*s2_11.*conj(H12) + H12.*s2_22.*conj(H22) + H13.*s2_33.*conj(H32);
S_CR13_conj = H11.*s2_11.*conj(H13) + H12.*s2_22.*conj(H23) + H13.*s2_33.*conj(H33);
S_CR23_conj = H21.*s2_11.*conj(H13) + H22.*s2_22.*conj(H23) + H23.*s2_33.*conj(H33);
% in case the cross-spectra do not work, wrap the conjugate products with
% an abs(); also, check the non diagonal results of Hhermit*SIGMA*Htrasp
% (symbolically)


% Compute Coherence
% Max ha detto che queste formule vanno bene, per adesso.
COH12 = sqrt(S_CR12.^2./(S_RR.*S_BP));
COH13 = sqrt(S_CR13.^2./(S_RR.*S_PAT));
COH23 = sqrt(S_CR23.^2./(S_BP.*S_PAT));

% TODO: Consider introducing Gain and Phases, if a proper trivariate form is
% found (otherwise, same getaround of coherences)

% Frequency Grid
r_LF = lf_hf(1) < FS & FS <= lf_hf(2); % % binary vector, true at LF indexes
r_HF = lf_hf(3) < FS & FS <= lf_hf(4); % % binary vector, true at HF indexes
r_VLF = 0.003 < FS & FS <= lf_hf(1);   % % binary vector, true at VLF indexes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   QUA C'ERANO TUTTI GLI INDICI INTEGRATI                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for f = length(fsamp):-1:1
            df=(fsamp(f)/2)/prec;
%             if ~isnan(fsamp(f))
                S_RR_TOT(f) = sum(S_RR(:,f).*df,'omitnan');
                S_RR_VLF(f) = sum(S_RR(r_VLF(:,f),f).*df,'omitnan');
                S_RR_LF(f) = sum(S_RR(r_LF(:,f),f).*df,'omitnan');
                S_RR_HF(f) = sum(S_RR(r_HF(:,f),f).*df,'omitnan');
                S_BP_TOT(f) = sum(S_BP(:,f).*df,'omitnan');
                S_BP_VLF(f) = sum(S_BP(r_VLF(:,f),f).*df,'omitnan');
                S_BP_LF(f) = sum(S_BP(r_LF(:,f),f).*df,'omitnan');
                S_BP_HF(f) = sum(S_BP(r_HF(:,f),f).*df,'omitnan');
                S_PAT_TOT(f) = sum(S_PAT(:,f).*df,'omitnan');
                S_PAT_VLF(f) = sum(S_PAT(r_VLF(:,f),f).*df,'omitnan');
                S_PAT_LF(f) = sum(S_PAT(r_LF(:,f),f).*df,'omitnan');
                S_PAT_HF(f) = sum(S_PAT(r_HF(:,f),f).*df,'omitnan');
                S_CR12_TOT(f) = sum(S_CR12(:,f).*df,'omitnan');
                S_CR12_VLF(f) = sum(S_CR12(r_VLF(:,f),f).*df,'omitnan');
                S_CR12_LF(f) = sum(S_CR12(r_LF(:,f),f).*df,'omitnan');
                S_CR12_HF(f) = sum(S_CR12(r_HF(:,f),f).*df,'omitnan');
                S_CR13_TOT(f) = sum(S_CR13(:,f).*df,'omitnan');
                S_CR13_VLF(f) = sum(S_CR13(r_VLF(:,f),f).*df,'omitnan');
                S_CR13_LF(f) = sum(S_CR13(r_LF(:,f),f).*df,'omitnan');
                S_CR13_HF(f) = sum(S_CR13(r_HF(:,f),f).*df,'omitnan');
                S_CR23_TOT(f) = sum(S_CR23(:,f).*df,'omitnan');
                S_CR23_VLF(f) = sum(S_CR23(r_VLF(:,f),f).*df,'omitnan');
                S_CR23_LF(f) = sum(S_CR23(r_LF(:,f),f).*df,'omitnan');
                S_CR23_HF(f) = sum(S_CR23(r_HF(:,f),f).*df,'omitnan');
                COH12_TOT(f) = sum(COH12(:,f).*df,'omitnan');
                COH12_VLF(f) = sum(COH12(r_VLF(:,f),f).*df,'omitnan');
                COH12_LF(f) = sum(COH12(r_LF(:,f),f).*df,'omitnan');
                COH12_HF(f) = sum(COH12(r_HF(:,f),f).*df,'omitnan');
                COH13_TOT(f) = sum(COH13(:,f).*df,'omitnan');
                COH13_VLF(f) = sum(COH13(r_VLF(:,f),f).*df,'omitnan');
                COH13_LF(f) = sum(COH13(r_LF(:,f),f).*df,'omitnan');
                COH13_HF(f) = sum(COH13(r_HF(:,f),f).*df,'omitnan');
                COH23_TOT(f) = sum(COH23(:,f).*df,'omitnan');
                COH23_VLF(f) = sum(COH23(r_VLF(:,f),f).*df,'omitnan');
                COH23_LF(f) = sum(COH23(r_LF(:,f),f).*df,'omitnan');
                COH23_HF(f) = sum(COH23(r_HF(:,f),f).*df,'omitnan');
%             end
end

% outs = {'GAIN_21_LF','GAIN_21_HF','PHASE_21_LF','PHASE_21_HF','GAIN_12_LF','GAIN_12_HF',...
%     'PHASE_12_LF','PHASE_12_HF','GCI_21_LF','GCI_21_HF','GCI_12_LF','GCI_12_HF',...
%     'COH_LF','COH_HF','F_LF','F_HF','S_RR_TOT','S_RR_VLF','S_RR_LF','S_RR_HF','S_BP_TOT','S_BP_VLF','S_BP_LF',...
%     'S_BP_HF','S_CR_TOT','S_CR_VLF','S_CR_LF','S_CR_HF',...
%     'S_RR','S_BP','S_CR','COH','GAIN_2to1','PHASE_2to1',...
%     'GAIN_1to2','PHASE_1to2','FS','r_LF','r_HF','GCI_2to1','GCI_1to2','H11','H12','H21','H22','S_CR_conj'};
  outs = {'FS','r_LF','r_HF','r_VLF','H11','H12','H13','H21','H22','H23','H31','H32','H33',...
          'S_RR','S_BP','S_PAT','S_CR12','S_CR13','S_CR23','S_CR12_conj','S_CR13_conj','S_CR23_conj',...
          'COH12','COH13','COH23',...
          'S_RR_TOT','S_RR_VLF','S_RR_LF','S_RR_HF',...
          'S_BP_TOT','S_BP_VLF','S_BP_LF','S_BP_HF',...
          'S_PAT_TOT','S_PAT_VLF','S_PAT_LF','S_PAT_HF',...
          'S_CR12_TOT','S_CR12_VLF','S_CR12_LF','S_CR12_HF',...
          'S_CR13_TOT','S_CR13_VLF','S_CR13_LF','S_CR13_HF',...
          'S_CR23_TOT','S_CR23_VLF','S_CR23_LF','S_CR23_HF',...
          'COH12_TOT','COH12_VLF','COH12_LF','COH12_HF',...
          'COH13_TOT','COH13_VLF','COH13_LF','COH13_HF',...
          'COH23_TOT','COH23_VLF','COH23_LF','COH23_HF'};

% for out=1:length(ret)
%     vout{out}=eval(outs{ret(out)});
% end

vout = {FS,r_LF,r_HF,r_VLF,H11,H12,H13,H21,H22,H23,H31,H32,H33,...
          S_RR,S_BP,S_PAT,S_CR12,S_CR13,S_CR23,S_CR12_conj,S_CR13_conj,S_CR23_conj,...
          COH12,COH13,COH23,...
          S_RR_TOT,S_RR_VLF,S_RR_LF,S_RR_HF,...
          S_BP_TOT,S_BP_VLF,S_BP_LF,S_BP_HF,...
          S_PAT_TOT,S_PAT_VLF,S_PAT_LF,S_PAT_HF,...
          S_CR12_TOT,S_CR12_VLF,S_CR12_LF,S_CR12_HF,...
          S_CR13_TOT,S_CR13_VLF,S_CR13_LF,S_CR13_HF,...
          S_CR23_TOT,S_CR23_VLF,S_CR23_LF,S_CR23_HF,...
          COH12_TOT,COH12_VLF,COH12_LF,COH12_HF,...
          COH13_TOT,COH13_VLF,COH13_LF,COH13_HF,...
          COH23_TOT,COH23_VLF,COH23_LF,COH23_HF};

% Qua c'era una serie di assegnamenti commentati, riferirsi a
% 'TriBaro_Redux.m'