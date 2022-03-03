function [Thetap1,Thetap2,Kappa1,Var2,opt] = regr_likel_cov(EKGR, COV, varargin)
% function [Thetap,Kappa,opt] = regr_likel(EKGR, COV, varargin)
% EKGR: R peaks times
% COV: Mx1 cell vector contaning Nx2 arrays where COV{1}(:,1) are times and
% COV {2}(:,2) are corresponding values.
% Thetap1: AR coefficient of IG AR model [(P*N)x1]
% Thetap2: AR coefficient of Gaussian AR model [(P*N)x(M+1)], 
% Thetap2(:,1) is for Gaussian model of RR, Thetap2(:,2:end) is for
% Gaussian model of COVs
%
% Copyright (C) Maximiliano Mollura, Luca Citi and Riccardo Barbieri, 2019-2020.
% All Rights Reserved. See LICENSE.TXT for license details.
% maximiliano.mollura@polimi.it
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
% uk = observ_ev(P+2:end);

% ck_t =[];
ck_v =[];
for i=1:length(COV)
%     tmp = COV{i}(1:end-1,1);
%     ck_t = [ck_t,tmp(:)];
    % Remove average value
    tmp = COV{i}(1:end,2);%-mean(COV{i}(1:end,2)); 
    ck_v = [ck_v,tmp(:)];
end

rr = diff(observ_ev);

wn = rr(P+1:end);
wcn = ck_v(P+1:end,:);

xn = toeplitz(rr(P:end-1), rr(P:-1:1));
cn = toeplitz(ck_v(P:end-1,1), ck_v(P:-1:1,1));
% Add others, if any
for i=2:length(COV)
    cn = [cn,toeplitz(ck_v(P:end-1,i), ck_v(P:-1:1,i))];
end

if opt.hasTheta0
    xn = [ones(length(wn),1), xn,cn];
end
[Thetap1,Kappa1,steps1,loglikel1] = maximize_loglikel(xn, wn);
[Thetap2,Var2] = WLS(xn, [wn,wcn]);
% Thetap2 = (xn'*xn)\xn'*[wn,wcn]; % Covariates Value (e.g. Systolic values) are Gaussian
% Var2 = (1/(size(xn,1)-P-1)).*sum(power([wn,wcn]-xn*Thetap2,2));

% Predict
opt.pred.Y_real = [wn,wcn];
opt.pred.Y_hat_IG = xn*Thetap1;
opt.pred.Y_hat_G = xn*Thetap2;

if opt.hasTheta0
    opt.Theta1_0 = Thetap1(1,:);
    Thetap1(1,:) = [];
    opt.Theta2_0 = Thetap2(1,:);
    Thetap2(1,:) = [];

end

opt.steps1 = steps1;
opt.loglikel1 = loglikel1;

% lsq_weighted