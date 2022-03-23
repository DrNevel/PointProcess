
function [Thetap,Mu,Kappa,L,opt] = pplikel_invgauss(EKGR, varargin)
% function [Thetap,Mu,Kappa,L,opt] = pplikel(EKGR, varargin)
%
%
% Copyright (C) Luca Citi and Riccardo Barbieri, 2010-2011.
% All Rights Reserved. See LICENSE.TXT for license details.
% {lciti,barbieri}@neurostat.mit.edu
% http://users.neurostat.mit.edu/barbieri/pphrv

% Default options
opt.delta = .005; % time increment in updating parameters (in seconds)
opt.P = 9; % RR order
opt.hasTheta0 = 1; % wether or not the AR model has a theta0 constant to account for the average mu
opt.weight = 0.98; % weighting factor
opt.W = 60; % window length for local likelihood estimate (in seconds)

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
delta = opt.delta;
W = opt.W;
P = opt.P;

opt.t0 = EKGR(1);
EKGR = EKGR(:) - opt.t0; % times are relative to EKGR(1) which is set to 0
J = floor(EKGR(end) / delta) + 1;

lastRi = find(EKGR(1:min(end,floor(5*W))) > W, 1) - 1; % index of last R event within W (check only first 5*W beats)
observ_ev = EKGR(1:lastRi);

Thetap = NaN(P + opt.hasTheta0, J);
Mu = NaN(1, J);
Kappa = NaN(1, J);
steps = zeros(1, J);
meanRR = NaN(1, J);
L = NaN(1, J);
opt.Loglikel = NaN(1, J);
flagnan = zeros(1,J);
k = 1000;

if isfield(opt,'thetap')
    thetap = opt.thetap;
else
    thetap = [];
end
fprintf('Finita parte matlab\n')

for j = ceil(W / delta):J
%     if mod(j,5000)==0
%         fprintf('iter %d \n',j);
%     end
    
    time = (j-1) * delta;
    if ~isempty(observ_ev) && (observ_ev(1) < time - W) 
        observ_ev(1) = []; % remove older event (there could be only one because we assume that in any delta interval there is at most one event)
        thetap = []; % force re-evaluation of starting point for thetap
    end
    event = EKGR(lastRi + 1) <= time; % whether an event happened in ((j-1)*delta,j*delta)
    if event
        lastRi = lastRi + 1;
        R = EKGR(lastRi);
        observ_ev(end+1,1) = R; % append current event
        thetap = []; % force re-evaluation of starting point for thetap
    end
    % if thetap is empty (i.e., observ_ev has changed) re-evaluate the variables that depend on observ_ev
    if isempty(thetap)
        uk = observ_ev(P+2:end); % length observ_ev - (P+2) 
        rr = diff(observ_ev); % length onbserv_ev-1
        wn = rr(P+1:end); % length rr - (P+1)
        xn = []; xt = [];
        if opt.hasTheta0
            xn = ones(length(wn),1);
            xt = 1;
        end
        xn = [xn, toeplitz(rr(P:end-1), rr(P:-1:1))];
        xt = [xt; rr(end:-1:end-P+1)];
        eta = weights(time, uk, opt.weight);
        %fprintf('ciclo thetap is empty\n')
        disp(j)

        % first version: funziona se uncensored (no thetap negli input)
        [thetap, k, loglikel] = maxi_loglike_invgauss(xn, wn, eta); % the uncensored loglikelihood is a good starting point
    else
        eta = weights(time, uk, opt.weight);
        
    end
    wt = time - observ_ev(end);
    
%     [L(j)] = Lambda_eval(xt, wt, thetap', k);
     
%      thetaTemp = thetap;
%      kTemp = k;
     
     %parte di right censoring: non funziona per ora
     %[thetap, k, loglikel, L] = maxi_loglikeRC(xn, wn, eta, thetap, k, xt, wt);
     
%     if(sum(isnan(thetap)) > 0 || sum(isnan(k)) > 0)
%         %fprintf('OH NO %d!\n',j);
%         thetap = thetaTemp;
%         k = kTemp;
%         flagnan(j) = 1;
%     end
    
    %disp(loglikel)
%     disp(thetap)
%     disp(k)
    
    %steps(j) = steps(j) + stepsj;
    mu = thetap * xt;
    Mu(j) = mu;
%     disp(mu);
    Thetap(:,j) = thetap;
    Kappa(j) = k;
    meanRR(j) = eta' * wn / sum(eta);
    opt.WT(j)=wt;
    opt.LogLikel(:,j) = sum(loglikel);
    
end
fprintf('uscito dal for\n')
if opt.hasTheta0
    opt.Theta0 = Thetap(1,:);
    Thetap(1,:) = [];
end

% opt.steps = steps;
opt.meanRR = meanRR;
opt.flagnan = flagnan;

