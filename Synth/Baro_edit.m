function [outs,vout] = Baro_edit(RR_Theta,COV_Theta, ord, s2_11, s2_22, fsamp, prec,ret,lf_hf)
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

if ~exist('lf_hf','var')
    lf_hf = [.04 .15; .15 .45]; %%%LF and HF ranges
elseif isempty(lf_hf)
    lf_hf = [.04 .15; .15 .45];
end

% Checks
[lt, N] = size(RR_Theta);
[lc, Nc] = size(COV_Theta);
assert(N == Nc);
assert(lt == lc);
assert(lt == 2*ord);


% Freqeuncy Axis
FS = linspace(0, 1/2, prec); % normalized fsamp = from 0:1/2
FS = FS'.*fsamp(:)';

% variances - cross-variances assumed to be null

% ARs Coefficient
A11 = fft(RR_Theta(1:ord,:),prec*2,1); %%% 1 --> return fft by columns
A12 = fft(RR_Theta(ord+1:end,:),prec*2,1);
A21 = fft(COV_Theta(1:ord,:),prec*2,1);
A22 = fft(COV_Theta(ord+1:end,:),prec*2,1);

A11 = A11(1:prec,:);
A12 = A12(1:prec,:);
A21 = A21(1:prec,:);
A22 = A22(1:prec,:);

H11 = (1-A22)./( ((1-A11).*(1-A22))-(A21.*A12) );
H22 = (1-A11)./( ((1-A11).*(1-A22))-(A21.*A12) );
H12 = (A12)  ./( ((1-A11).*(1-A22))-(A21.*A12) );
H21 = (A21)  ./( ((1-A11).*(1-A22))-(A21.*A12) );

% Compute Spectra
S_RR = abs(H11).^2 .* s2_11 + abs(H12).^2 .* s2_22;
S_BP = abs(H21).^2 .* s2_11 + abs(H22).^2 .* s2_22;
S_CR = abs(conj(H11).*H21) .* s2_11 + abs(conj(H12).*H22) .* s2_22;
S_CR_conj = abs(conj(H21).*H11) .* s2_11 + abs(conj(H22).*H12) .* s2_22;

outs = {'S_RR','S_BP','S_CR','S_CR','S_CR_conj','H11','H22','H12','H21'}; %ret =9

for out=1:length(ret)
    vout{out}=eval(outs{ret(out)});
end

end
