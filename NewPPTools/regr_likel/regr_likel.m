function [Thetap,Kappa,opt] = regr_likel(EKGR, varargin)
% function [Thetap,Kappa,opt] = regr_likel(EKGR, varargin)
%
%
% Copyright (C) Luca Citi and Riccardo Barbieri, 2010-2011.
% All Rights Reserved. See LICENSE.TXT for license details.
% {lciti,barbieri}@neurostat.mit.edu
% http://users.neurostat.mit.edu/barbieri/pphrv

% Default options
opt.P = 9; % RR order
opt.hasTheta0 = 1; % wether or not the AR model has a theta0 constant to account for the average mu
opt.maximize_loglikel = @likel_invnorm_mex; % use loglikelihood of inverse gaussian

% PROCESS OPTIONS
opt_names = fieldnames(opt);
J = 1;
while J <= length(varargin)
    key = varargin{J};
    keyN = find(strcmpi(key, opt_names));
    if isempty(keyN)
        warning('PPLIKEL:UnknownOption', 'Do not know option ''%s''', key);
    else
        opt.(opt_names{keyN}) = varargin{J+1};
    end
    J = J + 2;
end

% proxies for options
P = opt.P;
maximize_loglikel = opt.maximize_loglikel;

observ_ev = EKGR(:) - EKGR(1);
uk = observ_ev(P+2:end);
rr = diff(observ_ev);
wn = rr(P+1:end);
xn = toeplitz(rr(P:end-1), rr(P:-1:1));
if opt.hasTheta0
    xn = [ones(length(wn),1), xn];
end
[Thetap,Kappa,steps,loglikel] = maximize_loglikel(xn, wn);
[Thetap2,Var2] = WLS(xn, wn);
Var2 = sum(Var2);
% Thetap2 = (xn'*xn)\xn'*wn; % Least Square Solution
% Var2 = (1/(size(xn,1)-P-1)).*sum(power(wn-xn*Thetap2,2));

% Predict
opt.pred.Y_real = wn;
opt.pred.Y_hat_IG = xn*Thetap;
opt.pred.Y_hat_G = xn*Thetap2;

if opt.hasTheta0
    opt.Theta0 = Thetap(1,:);
    Thetap(1,:) = [];
    opt.Gauss.Theta0 = Thetap2(1,:);
    Thetap2(1,:) = [];
end

opt.steps = steps;
opt.loglikel = loglikel;
opt.Gauss.Thetap = Thetap2;
opt.Gauss.Var = Var2;

