function [Thetap,Mu,Kappa,L,opt] = pplikel_cov_corr_multi(EKGR, COV, varargin)
% function [Thetap,Mu,Kappa,L,opt] = pplikel_cov(EKGR, COV, varargin)
% Example: [Thetapc,Muc,Kappac,Lc,optc] = pplikel_cov(EKGR, {[COV_T(:), COV_VAL(:)]},...
%             'P', 4, 'hasTheta0', 1, 'delta', 0.005,'W', 60, 'weight', 0.98);
% Example: [Thetapc,Muc,Kappac,Lc,optc] = pplikel_cov(EKGR, {[COV_T(:), COV_VAL(:)]},...
%             'P', 4, 'hasTheta0', 1, 'delta', 0.005,'W', 60, 'weight', 0.98,'P2',9);
% EKGR: R-Peaks, Nx1 elements
% COV: cell array Mx1, each with Nx2 elements, times and values
% P: AR order, final order will be P*(M+1) (def: 4)
% hasTheta0: flag to include a coeff for 1 (def: 1)
% delta: time resolution (def: 0.005)
% W: prediction window (def:60)
% weight: weighting factor (def:0.98)
% P2: order of monovariate AR models on COV only (def:[])
%
% Copyright (C) Maximiliano Mollura, Luca Citi and Riccardo Barbieri, 2019-2020.
% All Rights Reserved. See LICENSE.TXT for license details.
% maximiliano.mollura@polimi.it
% riccardo.barbieri@polimi.it
% http://users.neurostat.mit.edu/barbieri/pphrv

% Default options
opt.delta = .005; % time increment in updating parameters (in seconds)
opt.P = 4; % RR-COV order
opt.P2 = []; % COV only order
opt.hasTheta0 = 1; % wether or not the AR model has a theta0 constant to account for the average mu
opt.weight = 0.98; % weighting factor
opt.W = 60; % window length for local likelihood estimate (in seconds)
opt.maximize_loglikel = @likel_invnorm_mex; % use loglikelihood of inverse gaussian
opt.Rsim = [];
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
P2 = opt.P2;
if ~isempty(P2)
    mono_cov_flg = 1;
else
    mono_cov_flg = 0;
end

maximize_loglikel = opt.maximize_loglikel;

opt.t0 = EKGR(1);
EKGR = EKGR(:) - opt.t0; % times are relative to EKGR(1) which is set to 0
J = floor(EKGR(end) / delta) + 1;

lastRi = find(EKGR(1:min(end,floor(5*W))) > W, 1) - 1; % index of last R event within W (check only first 5*W beats)

if isempty(opt.Rsim)
    observ_ev = EKGR(1:lastRi);
    uk0 = NaN;
else
    observ_ev = opt.Rsim(1:lastRi); % To run simulations where EKGR is a time at a fixed freq.
    uk0 = EKGR(1:lastRi);
end

S=size(COV{1});
ck_t = NaN(S(2)-1,S(1));
ck_v = NaN(S(2)-1,S(1));
for i=length(COV):-1:1 % COV{i} is assumed to have same length of EKGR
    tmp_cov = COV{i}';
    tmp = tmp_cov(1:end-1,1) - opt.t0;
    ck_t(:,i) = tmp(:);
    % Remove average value
    tmp = tmp_cov(1:end-1,2);%-mean(COV{i}(1:end-1,2)); % The last Cov value influences the successive RR, so we remove it.
    ck_v(:,i) = tmp(:);
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
meanCOV = NaN(length(COV), J);
opt.LogLikel = NaN(1, J);
Thetap2 = NaN(P*(length(COV)+1) + opt.hasTheta0, length(COV)+1, J);
Var2 =  NaN(length(COV)+1, J);
Thetap_Cov = NaN(P2 + opt.hasTheta0, length(COV), J);
Mu_Gaus = NaN(1, J);
Cov_Gaus = NaN(length(COV), J);
if mono_cov_flg
    Mu_Cov = NaN(length(COV), J);
    Var_Cov =  NaN(length(COV), J);
end
thetap = [];

time_ax = 0:delta:EKGR(end);

for j = ceil(W / delta):J
    
    time = (j-1) * delta;
    if ~isempty(observ_ev) && (observ_ev(1) < time - W) && isempty(opt.Rsim)
        observ_ev(1) = []; % remove older event (there could be only one because we assume that in any delta interval there is at most one event)
        thetap = []; % force re-evaluation of starting point for thetap
        cc(1,:) = [];
        cc_t(1,:) = [];
    elseif ~isempty(observ_ev) && (uk0(1) < time - W) && ~isempty(opt.Rsim)
        observ_ev(1) = []; % remove older event (there could be only one because we assume that in any delta interval there is at most one event)
        uk0(1) = [];
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
            uk0(end+1) = EKGR(lastRi);
        end
        observ_ev(end+1,1) = R; % append current event
        C = ck_v(lastRi-1,:);
        CT = ck_t(lastRi-1,:);
        cc(end+1,:) = C; % Ideally append it when encounter not at R
        cc_t(end+1,:) = CT;
        thetap = []; % force re-evaluation of starting point for thetap
    end
    % if thetap is empty (i.e., observ_ev has changed) re-evaluate the variables that depend on observ_ev
    if isempty(thetap)
        
        if isempty(opt.Rsim)
            uk = observ_ev(P+2:end);
            rr = diff(observ_ev);
        else
            uk = uk0(P+1:end-1);
            rr = observ_ev(1:end-1); % To run simulations where EKGR is a time at a fixed freq.
        end
        wn = rr(P+1:end); % RR outcomes
        wcn = cc(P+1:end,:); % COV outcomes

        for i=length(COV):-1:1 % COV{i} is assumed to have same length of EKGR
            clear wn2
            % Extract Positions for Mu regression
            [~,min_ix] = min(abs(time_ax-cc_t(1,i)));
            [~,max_ix] = min(abs(time_ax-cc_t(end,i)));
            [~,mu_cc_t(:,i)] = min(abs(time_ax(min_ix:max_ix)-cc_t(:,i)),[],2);
            wn2(:,i) = Mu(mu_cc_t(:,i)+min_ix);
            clear mu_cc_t
        end
        wn2 = wn2(P+1:end,:);
        
        cn = toeplitz(cc(P:end-1,1), cc(P:-1:1,1)); % Stop at end-1 to use the last to train the censoring part
        ct = cc(end:-1:end-P+1,1);
        
        for cl = 2:length(COV) % values always after times
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
        %         [thetap2,var2] = WLS(xn, [wn,wcn], eta);
        
    else
        eta = weights(time, uk, opt.weight);
    end
    
    wt = time - observ_ev(end); % THESE ARE THE TIMES OF OUTCOMES FOR CENSORING!
    [thetap, k, stepsj, L(j), loglikel] = maximize_loglikel(xn, wn, eta, thetap, k, xt, wt);
    if sum(isnan(wn2))~=0
        [thetap2,var2] = WLS(xn, [wn,wcn], eta); % here the weights 'eta' are changing
    else
        [thetap2,var2] = WLS(xn, [wn2,wcn], eta); % here the weights 'eta' are changing
    end

    
    if mono_cov_flg
        for cl = 1:length(COV)
            % Covarite AR on COV only
            wcn2 = cc(P2+1:end,cl); % COV outcomes
            cn2 = toeplitz(cc(P2:end-1,cl), cc(P2:-1:1,cl)); % Stop at end-1 to use the last to train the censoring part
            eta2 = weights(time, observ_ev(P2+2:end), opt.weight); % Probably better to use the pressure times but here using R times
            ct2(:,cl) = [1;cc(end:-1:end-P2+1,cl)];
            [thetap_cov(:,cl),cov_var(:,cl)] = WLS([ones(size(cn2,1),1),cn2], wcn2, eta2);
        end
    end
    
    if (sum(isnan(thetap))>=1 || sum(isnan(k)) > 0)
        fprintf('\n NaN coeffs at time: %d [sec]',round(time,2))
        if j>1 && sum(isnan(Thetap(:,j-1)))<1
            [thetap, k, steps(j)] = maximize_loglikel(xn, wn, eta); % the uncensored loglikelihood is a good starting point
            if (sum(isnan(thetap))>=1 || sum(isnan(k)) > 0)
                thetap = Thetap(:,j-1);
            end
        end
    end
    steps(j) = steps(j) + stepsj;
    mu = thetap' * xt;
    
    % IG Params
    Mu(j) = mu;
    Thetap(:,j) = thetap;
    Kappa(j) = k;
    meanRR(j) = eta' * wn / sum(eta);
    tmp = eta' * wcn / sum(eta);
    meanCOV(:,j) = tmp';
    opt.LogLikel(:,j) = sum(loglikel);
    % Gauss Params
    Mu_Gaus(j) = thetap2(:,1)'*xt; % RR Bivariate series Gaussian
    Cov_Gaus(:,j) = thetap2(:,2:end)'*xt; % COV Bivariate series Gaussian
    Thetap2(:,:,j) = thetap2;
    Var2(:,j) = var2';
    %     if mono_cov_flg
    %         % Cov Params
    %         Mu_Cov(:,j) = diag(thetap_cov'*ct2); % COV Monovariate series Gaussian
    %         Thetap_Cov(:,:,j) = thetap_cov;
    %         Var_Cov(:,j) = cov_var;
    %     end
    
end

if opt.hasTheta0
    opt.Theta0 = Thetap(1,:);
    opt.Theta2_0 = permute(Thetap2(1,:,:),[2,3,1]);
    Thetap(1,:) = [];
    Thetap2(1,:,:) = [];
    if ~isempty(P2)
        Thetap_Cov(1,:,:) = [];
    end
end

opt.steps = steps;
opt.meanRR = meanRR;
opt.meanCOV = meanCOV;

opt.Mu_Gaus = Mu_Gaus;
opt.Cov_Gaus = Cov_Gaus;
opt.Thetap2 = permute(Thetap2,[1,3,2]);
opt.Var2 = Var2;
if mono_cov_flg
    opt.Mu_Cov = Mu_Cov;
    opt.Thetap_cov = permute(Thetap_Cov,[1,3,2]);
    opt.Var_Cov = Var_Cov;
end