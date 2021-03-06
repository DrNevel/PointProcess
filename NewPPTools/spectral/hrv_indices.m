function [powLF, powHF, bal, warn, powVLF, powTot, ToT, F] = hrv_indices(Thetap, Var, fsamp, prec, dwsample)
% function [powLF, powHF, bal, warn] = hrv_indices(Thetap, Var, fsamp, dwsample)
%
% Evaluate Heart Rate Variability indices from the outputs of
% the point-process likelihood modeling (pplikel.m)
%
%
% Copyright (C) Luca Citi and Riccardo Barbieri, 2010-2011.
% All Rights Reserved. See LICENSE.TXT for license details.
% {lciti,barbieri}@neurostat.mit.edu
% http://users.neurostat.mit.edu/barbieri/pphrv

[P,J] = size(Thetap);

% initialize output
% powLF and powHF are NaN where they cannot be evaluated (e.g. at the beginning)
powLF = NaN(1,J);
powHF = NaN(1,J);
powVLF = NaN(1,J);
powTot = NaN(1,J);
if exist('prec','var')
    ToT = NaN(prec,J);
    F = NaN(prec,J);
else
    prec = 512;
end
% warn==1 is OK
% warn~=1 warning:
%   mod(warn, 2) is the scale factor used to shrink the poles in case of instability (slightly less than 1 is fine)
%   bitand(floor(warn), 2) if powLF was negative (increase AR order?)
%   bitand(floor(warn), 4) if powHF was negative (increase AR order?)
warn = NaN(1,J);

for i = 1:J % For each time instant
    if isnan(Thetap(1,i))
        continue;
    end
    [tot,comps,compsp,f,pole_freq,pole_pow,pole_res,poles,mod_scale] = spectral(Thetap(:,i), Var(i), fsamp(i), prec, 0);
    warn(i) = mod_scale;
    pf = abs(pole_freq);
    powLF(i) = sum(pole_pow((pf>0.04) & (pf<0.15)));
    if powLF(i) <= 0
        powLF(i) = powLF(i-1);
        warn(i) = warn(i) + 2;
    end
    powHF(i) = sum(pole_pow((pf>0.15) & (pf<0.45)));
    if powHF(i) <= 0
        powHF(i) = powHF(i-1);
        warn(i) = warn(i) + 4;
    end
    powVLF(i) = sum(pole_pow(pf<0.04));
    powTot(i) = sum(pole_pow);
    
    if nargout>6
        ToT(:,i) = tot;
        F(:,i) = f;
    end
end

b = ceil(J/10);
if b > 1
    b = hamming(min(21, b));
    b = b/sum(b);
    powLF = fliplr(filter(b,1,fliplr(filter(b,1,powLF))));
    powHF = fliplr(filter(b,1,fliplr(filter(b,1,powHF))));
    powVLF = fliplr(filter(b,1,fliplr(filter(b,1,powVLF))));
    powTot = fliplr(filter(b,1,fliplr(filter(b,1,powTot))));
end

if nargin > 4 && dwsample ~= 1
    powLF = decimate(powLF, dwsample, 'fir');
    powHF = decimate(powHF, dwsample, 'fir');
    powVLF = decimate(powVLF, dwsample, 'fir');
    powTot = decimate(powTot, dwsample, 'fir');
end

bal = powLF ./ powHF;
