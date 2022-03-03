function [EKGR,SYST,DIAST,ONSET,RR_correction] = closeECG(EKGR,SYST,DIAST,ONSET,correct_diast)

% Copyright (C) Maximiliano Mollura and Riccardo Barbieri, 2019-2020.
% All Rights Reserved. See LICENSE.TXT for license details.
% maximiliano.mollura@polimi.it
% barbieri@neurostat.mit.edu

if ~exist('SYST','var')
elseif isempty(SYST),clear SYST,end
if ~exist('DIAST','var')
elseif isempty(DIAST),clear DIAST,end
if ~exist('ONSET','var')
elseif isempty(ONSET),clear ONSET,end

hrv = diff(EKGR);
correctECG = find(hrv>=3|hrv<0.3);
RR_val = hrv(correctECG);
RR_correction = [correctECG(:), RR_val(:), RR_val(:)];
mdn = NaN(length(correctECG),1);
k = 1;
% fare correzione pressione
% correctECG = sort(correctECG,'descend');
for i = correctECG
    if (i - 5 > 0) && (i + 5 <= length(hrv))
        mdn(k) = median(hrv(i-5:i+5));
    elseif (i - 5 <= 0) && (i + 10 <= length(hrv))
        mdn(k) = median(hrv(i:i+10));
    elseif (i + 5 > length(hrv)) && (i - 10 > 0)
        mdn(k) = median(hrv(i-10:i));
    else
        mdn(k) = median(hrv(i-10:i));
    end
    EKGR(i+1:end) = EKGR(i+1:end)-(EKGR(i+1)-EKGR(i))+ mdn(k);
    if exist('SYST','var')
        SYST(i+1:end) = SYST(i+1:end)-(SYST(i+1)-SYST(i))+ mdn(k);
    end
    if exist('DIAST','var') && correct_diast==0 %%%
        DIAST(i+2:end) = DIAST(i+2:end)-(DIAST(i+2)-DIAST(i+1))+ mdn(k);
    elseif exist('DIAST','var') && correct_diast==1
        DIAST(i+1:end) = DIAST(i+1:end)-(DIAST(i+1)-DIAST(i))+ mdn(k);
    end
    if exist('ONSET','var')
        ONSET(i+1:end) = ONSET(i+1:end)-(ONSET(i+1)-ONSET(i))+ mdn(k);
    end
    k = k+1;
end
RR_correction(:,1) = EKGR(correctECG(:)); % switch to seconds
RR_correction(:,3) = mdn;

if ~exist('SYST','var'),SYST=[];end
if ~exist('DIAST','var'),DIAST=[];end
if ~exist('ONSET','var'),ONSET=[];end
