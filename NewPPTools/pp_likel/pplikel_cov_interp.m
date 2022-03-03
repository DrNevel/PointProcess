function [Thetap,Mu,Kappa,L,opt] = pplikel_cov(EKGR, COV, varargin)
% function [Thetap,Mu,Kappa,L,opt] = pplikel(EKGR, varargin)
%
%
% Copyright (C) Maximiliano Mollura, Luca Citi and Riccardo Barbieri, 2019-2020.
% All Rights Reserved. See LICENSE.TXT for license details.
% {lciti,barbieri}@neurostat.mit.edu
% http://users.neurostat.mit.edu/barbieri/pphrv

% Default options
opt.delta = .005; % time increment in updating parameters (in seconds)
opt.P = 9; % RR order
opt.hasTheta0 = 1; % wether or not the AR model has a theta0 constant to account for the average mu
opt.weight = 0.98; % weighting factor
opt.W = 60; % window length for local likelihood estimate (in seconds)
opt.maximize_loglikel = @likel_invnorm_mex; % use loglikelihood of inverse gaussian
opt.Rsim = [];
opt.f_res = [];

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
maximize_loglikel = opt.maximize_loglikel;

opt.t0 = EKGR(1);
EKGR = EKGR(:) - opt.t0; % times are relative to EKGR(1) which is set to 0
J = floor(EKGR(end) / delta) + 1;

lastRi = find(EKGR(1:min(end,floor(5*W))) > W, 1) - 1; % index of last R event within W (check only first 5*W beats)
if isempty(opt.Rsim)
    observ_ev = EKGR(1:lastRi);
else
    observ_ev = opt.Rsim(1:lastRi); % To run simulations where EKGR is a time at a fixed freq.
end
ck_t =[];
ck_v =[];
for i=1:length(COV) % COV{i} is assumed to have same length of EKGR
    tmp = COV{i}(1:end-1,1) - opt.t0;
    ck_t = [ck_t,tmp(:)];
    % Remove average value
    tmp = COV{i}(1:end-1,2)-mean(COV{i}(1:end-1,2)); % The last Cov value influences the successive RR, so we remove it.
    ck_v = [ck_v,tmp(:)];
end
cc = ck_v(1:lastRi-1,:);
cc_t = ck_t(1:lastRi-1,:);

% Init
Thetap = NaN(P*(length(COV)+1) + opt.hasTheta0, J);
Mu = NaN(1, J);
Kappa = NaN(1, J);
steps = zeros(1, J);
L = NaN(1, J);
meanRR = NaN(1, J);
opt.LogLikel = NaN(1, J);
Thetap2 = NaN(P*(length(COV)+1) + opt.hasTheta0, length(COV)+1, J);
Var2 =  NaN(length(COV)+1, J);
Mu2 = NaN(1, J);
COV2 = NaN(1, J);
thetap = [];

for j = ceil(W / delta):J

    time = (j-1) * delta;    
    if ~isempty(observ_ev) && (observ_ev(1) < time - W)
        observ_ev(1) = []; % remove older event (there could be only one because we assume that in any delta interval there is at most one event)
        thetap = []; % force re-evaluation of starting point for thetap
        cc(1,:) = [];
        cc_t(1,:) = [];
    end

    event = EKGR(lastRi + 1) <= time; % whether an event happened in ((j-1)*delta,j*delta]
    if event
        if mod(lastRi,100)==0
            fprintf('Processed Beats: %d/%d \n',lastRi,length(EKGR))
        end
        lastRi = lastRi + 1;
        if isempty(opt.Rsim)
            R = EKGR(lastRi);
        else
            R = opt.Rsim(lastRi); % To run simulations where EKGR is a time at a fixed freq.
        end
        observ_ev(end+1,1) = R; % append current event
        C = ck_v(lastRi-1,:);
        CT = ck_t(lastRi-1,:);
        cc(end+1,:) = C; % Ideally, append it when encounter, not at R
        cc_t(end+1,:) = CT;
        thetap = []; % force re-evaluation of starting point for thetap
    end
    % if thetap is empty (i.e., observ_ev has changed) re-evaluate the variables that depend on observ_ev
    if isempty(thetap)
        
        uk = observ_ev(P+2:end);        
        rr = diff(observ_ev);
        
        % Interpolation if requested within the window
        if ~isempty(opt.f_res)
            rr = pchip(observ_ev(2:end),rr,observ_ev(2):1/opt.f_res:observ_ev(end)); % Interpolate RR and then interpolate with this number of sample
            rr = rr(:);
            uk = linspace(uk(1),uk(end),length(rr)-P);
            uk = uk(:);
            for v = 1:size(cc,2)
                cc_int(:,v) = pchip(cc_t(:,v),cc(:,v),linspace(cc_t(1,v),cc_t(end,v),length(rr)));
                cc_t_int(:,v) = linspace(cc_t(1,v),cc_t(end,v),length(rr));
            end
            cc = cc_int;
            cc_t = cc_t_int;
            cc_int = [];
            cc_t_int = [];
        end
        
        wcn = cc(P+1:end,:); % COV outcomes
        wn = rr(P+1:end); % RR outcomes
        
        cn = toeplitz(cc(P:end-1,1), cc(P:-1:1,1)); % Stop at end-1 to use the last to train the censoring part
        ct = cc(end:-1:end-P+1,1);
        for cl = 2:length(COV)
            cn = [cn,toeplitz(cc(P:end-1,cl), cc(P:-1:1,cl))];
            ct = [ct;cc(end:-1:end-P+1,cl)];
        end
        
        xn = []; xt = [];
        if opt.hasTheta0
            xn = ones(length(wn),1);
            xt = 1;
        end

        xn = [xn, toeplitz(rr(P:end-1), rr(P:-1:1)),cn];
        xt = [xt; rr(end:-1:end-P+1);ct]; % for the censoring (it includes the last RR (rr(end))) THESE ARE OUTCOMES FOR CENSORING!
        eta = weights(time, uk, opt.weight);
        [thetap, k, steps(j)] = maximize_loglikel(xn, wn, eta); % the uncensored loglikelihood is a good starting point
        [thetap2,var2] = WLS(xn, [wn,wcn], eta);
    else
        eta = weights(time, uk, opt.weight);
    end
    wt = time - observ_ev(end); % THESE ARE THE TIMES OF OUTCOMES FOR CENSORING!
    [thetap, k, stepsj, L(j), loglikel] = maximize_loglikel(xn, wn, eta, thetap, k, xt, wt);
    [thetap2,var2] = WLS(xn, [wn,wcn], eta); % here the weights 'eta' are changing
    if sum(isnan(thetap))>=1
        fprintf('\n NaN coeffs at time: %d [sec]',round(time,2))
        if j>1 && sum(isnan(Thetap(:,j-1)))<1
            thetap = Thetap(:,j-1);
        end
    end
    steps(j) = steps(j) + stepsj;
    mu = thetap' * xt;

    % IG Params
    Mu(j) = mu;
    Thetap(:,j) = thetap;
    Kappa(j) = k;
    meanRR(j) = eta' * wn / sum(eta);
    opt.LogLikel(:,j) = sum(loglikel);
    % G Params
    Mu_Gaus(j) = thetap2(:,1)'*xt;
    Cov_Gaus(j) = thetap2(:,2)'*xt;
    Thetap2(:,:,j) = thetap2;
    Var2(:,j) = var2;

end

if opt.hasTheta0
    opt.Theta0 = Thetap(1,:);
    opt.Theta2_0 = Thetap2(1,:,:);
    Thetap(1,:) = [];
    Thetap2(1,:,:) = [];
end

opt.steps = steps;
opt.meanRR = meanRR;

opt.Mu2 = Mu2;
opt.COV2 = COV2;
opt.Thetap2 = permute(Thetap2,[1,3,2]);
opt.Var2 = Var2;