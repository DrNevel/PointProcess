function [outs_lbl,varargout] = Baro_vect(RR_Theta,COV_Theta, ord, s2, cov_s2, fsamp, prec, mode, ret,S_RR_mono,S_BP_mono)
%[lbl,GAIN_21_LF,GAIN_21_HF,PHASE_21_LF,PHASE_21_HF,GAIN_12_LF,GAIN_12_HF,...
%     PHASE_12_LF,PHASE_12_HF,GCI_21_LF,GCI_21_HF,GCI_12_LF,GCI_12_HF,...
%     COH_LF,COH_HF,F_LF,F_HF,S_RR_TOT,S_RR_VLF,S_RR_LF,S_RR_HF,S_BP_TOT,S_BP_VLF,S_BP_LF,S_BP_HF,S_CR_TOT,S_CR_VLF,S_CR_LF,S_CR_HF,...
%     S_RR,S_BP,S_CR,COH,GAIN_2to1,PHASE_2to1,...
%     GAIN_1to2,PHASE_1to2,FS,r_LF,r_HF,GCI_2to1,GCI_1to2,H11,H12,H21,H22,S_CR_conj] = Baro_vect(Feat.EKGR_1.PP.Cov.Thetap,Feat.EKGR_1.PP.Cov.covariates_ar_ff, PPord, Feat.EKGR_1.PP.Cov.s2,...
%     Feat.EKGR_1.PP.Cov.covariates_s2_ff, 1./Feat.EKGR_1.PP.Cov.meanRR, 512, 1, 1:29); %  mode set to 1 instead of 1:4
% 
% Last Update: 15/5/2020
%
% Copyright (C) Maximiliano Mollura and Riccardo Barbieri, 2019-2020.
% All Rights Reserved. See LICENSE.TXT for license details.
% maximiliano.mollura@polimi.it
% riccardo.barbieri@polimi.it

% Checks
[lt, N] = size(RR_Theta);
[lc, Nc] = size(COV_Theta);
assert(N == Nc);
assert(lt == lc);
assert(lt == 2*ord);

if isscalar(fsamp)
    fsamp = repmat(fsamp,[1,N]);
end

%%%ret = number of desired outputs
for i = 1:length(ret)
    varargout{i} = []; %1x46 cell
end

if N > 1024
    blk = 0:1024:N; % do it in blocks, 1024 is a reasonable vectorization/memory tradeoff
    %%% blk = series of blocks from 0 to N with size 1024
    blk(end) = N;
else
    blk = [0 N];
end

%%% length blk = ceil(N/1024)
%%% j from length blk to 2, backward
for j = length(blk):-1:2
    fprintf('\n Block %d / %d', length(blk)-j,length(blk))
    bb = blk(j-1)+1:blk(j);
    % bb = all unit steps inside each block; size 1024 except for the last block (first iteration)

    if exist('S_RR_mono','var')&&exist('S_BP_mono','var')
        [outs_lbl, out] = ...
            Baro(RR_Theta(:,bb),COV_Theta(:,bb), ord, s2(:,bb), cov_s2(:,bb), fsamp(:,bb), prec, mode,ret,[],S_RR_mono(:,bb),S_BP_mono(:,bb));
    elseif exist('S_RR_mono','var')&&~exist('S_BP_mono','var')
        [outs_lbl, out] = ...
            Baro(RR_Theta(:,bb),COV_Theta(:,bb), ord, s2(:,bb), cov_s2(:,bb), fsamp(:,bb), prec, mode,ret,[],S_RR_mono(:,bb));
        %RR_Theta(:,bb) extracts bb columns of RR_Theta
    else
        [outs_lbl, out] = ...
            Baro(RR_Theta(:,bb),COV_Theta(:,bb), ord, s2(:,bb), cov_s2(:,bb), fsamp(:,bb), prec, mode,ret);
            % s2_11 = s2
            % s2_22 = cov_s2 --> Features.PP.Cov.cov_s2(2,:)
    end
    for i = 1:length(ret)
        varargout{i}(:,bb) = out{i};
    end
end

outs_lbl = outs_lbl(ret);