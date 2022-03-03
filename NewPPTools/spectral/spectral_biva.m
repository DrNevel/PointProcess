function [tot,comps,compsp,f,pole_freq,pole_pow,pole_res,poles,mod_scale] = spectral_biva(A, Var, fsamp, f, do_comps,B,CovVar)
% function [tot,comps,compsp,f,pole_freq,pole_pow,pole_res,poles] = spectral(A, Var, fsamp, f, do_comps)
%
% Evaluate spectra and components
%
% [tot,comps,compsp,f,pole_freq,pole_pow,pole_res,poles] = spectral(A, Var, fsamp, f, do_comps)
% with:
%   A is a vector of length P with the regression terms of the IG model (Thetap)
%   Var is the variance of the IG
%   fsamp is the sampling frequency
%   f is the number of frequency points to use in [0 fsamp/2] if f positive or
%          [-fsamp/2 fsamp/2] if f negative
%   do_comps can be
%          0: do not find spectral components but only the total PSD
%          1: find spectral components both for each individual pole (output
%                    "comps") and for pair of conjugate poles (output "compsp")
%          2: as in the case 1 and also plot "compsp" (default)
%          3: as in the case 1 and also plot "comps"
% the output arguments are:
%   tot is a vector with the total PSD (plotted in black)
%   comps is a matrix with P rows, each row is the PSD of a single pole
%          (plotted in colour with do_comps=3)
%   compsp is a matrix with <=P rows, each row is the PSD of a single real pole
%          or of a pair of conjugate poles (plotted in colour with do_comps=2)
%   f is a vector with the frequency points where the PSD is evaluated
%   pole_freq is a vector with the nominal frequency of each pole
%   pole_pow is a vector with the power of each pole
%   pole_res is a vector with the residual of each pole
%   poles is a vector with the position of each pole in the complex plane
%
%
% Copyright (C) Maximiliano Mollura, Luca Citi and Riccardo Barbieri, 2019-2020.
% All Rights Reserved. See LICENSE.TXT for license details.
% maximiliano.mollura@polimi.it
% {lciti,barbieri}@neurostat.mit.edu
% http://users.neurostat.mit.edu/barbieri/pphrv


if nargin < 4
    f = 1024;
end
if nargin < 5
    do_comps = 2;
end


[P,N] = size(A);
assert(N==1)


% find poles
poles = roots([1; -A]);
[ix,ix] = sort(abs(angle(poles)));
poles = poles(ix);

if exist('B','var')
    if ~isempty(B)
        % find zeros of the covariate to correct the spectrum
        [P2,N] = size(B);
        assert(N==1)
        assert(P==P2)
        covs = roots([B;0]);
        [ix,ix] = sort(abs(angle(covs)));
        covs = covs(ix);
        cov_res = zeros(P2,1);
    end
end

% fix AR models that might have become slightly unstable due to the estimation process
% using an exponential decay (see Stoica and Moses, Signal Processing 26(1) 1992)
mod_scale = min(.99/max(abs(poles)), 1);
poles = poles * mod_scale;
A = A .* cumprod(ones(size(A))*mod_scale);

% find residuals
pole_res = zeros(P,1);
pole_freq = zeros(P,1);
pole_pow = zeros(P,1);
for i = 1:P
    rm = poles;
    rm(i) = [];
    ri = poles(i);
    pole_res(i) = 1 / (ri * prod(ri - rm) * prod(1/ri - conj(poles))); % Effect of the current pole + Residual effect of the other poles on the frequency of the current pole
        
    pole_freq(i, :) = angle(ri)/(2*pi)* fsamp;
%     pole_freq(i, :) = angle(ri)/(2*pi); % O riscali lo spazio Re-Im o
%     riscali i limiti degli indici di LF e HF.

     pole_pow(i, :) = Var * real(pole_res(i));

     if exist('B','var')
        if ~isempty(B)
            cov_res(i) = (prod(ri-covs)*prod(ri-conj(covs))) / (ri * prod(ri - rm) * prod(1/ri - conj(poles))); % Effect of the current pole + Residual effect of the other poles on the frequency of the current pole
            pole_pow(i, :) = Var * real(pole_res(i))+ CovVar * real(cov_res(i));
        end
     end
end

tot = [];
comps = [];
compsp = [];

if ~isempty(f)
    % components
    if length(f) > 1
        FS = f / fsamp; % normalized freq
    else
        FS = linspace((sign(f)-1)/4, 1/2, abs(f)); % normalized freq [-.5,.5] if f<0 else [0,.5]
%         FS = linspace((sign(f)-1)/4, fsamp/2, abs(f)); % normalized freq [-.5,.5] if f<0 else [0,.5]
%         FS = linspace(0, fsamp/2, abs(f)); % normalized freq [-.5,.5] if f<0 else [0,.5]
    end
%     z = exp((2j*pi*FS).*fsamp); % Moltiplicare per fsamp ci fa solo vedere la replica dello spettro, non aumenta la f_Nyq
%     z = exp((2j*pi*FS)).*fsamp; % NON HA SENSO
    z = exp(2j*pi*FS); % ORIGINALE
    if do_comps
        comps = zeros(P, length(z));
        refpoles = 1./conj(poles);
        for i = 1:P
            pp = pole_res(i) * poles(i) ./ (z - poles(i));
            refpp = -conj(pole_res(i)) * refpoles(i) ./ (z - refpoles(i));
            comps(i,:) = Var/fsamp * (pp + refpp);
        end
        compsp = comps(1, :);
        for i = 2:P
            if abs(poles(i) - conj(poles(i-1))) / abs(poles(i) + conj(poles(i-1)) + eps) < 1e-6
                compsp(end,:) = compsp(end,:) + comps(i, :);
            else
                compsp(end+1,:) = comps(i, :);
            end
        end
    end
    
    tot = Var/fsamp ./ abs(polyval([1; -A], conj(z))).^2; %% la frequenza di Nyquist è intrinsecamente definita all'interno dei coefficienti istante per istante
    f = FS * fsamp;
    
    if exist('B','var')
        if ~isempty(B)
                tot = Var/fsamp ./ abs(polyval([1; -A], conj(z))).^2 + CovVar/fsamp .* (abs(polyval([B;0], conj(z))).^2./abs(polyval([1; -A], conj(z))).^2); %% la frequenza di Nyquist è intrinsecamente definita all'interno dei coefficienti istante per istante
        end
    end

end
if do_comps >= 2
    figure;
    if do_comps == 2
        plot(f, real(compsp));
    else
        plot(f, real(comps));
    end
    hold on;
    plot(f, tot, 'k');
    xlim([f(1), f(end)]);
end

