function [EKGR,Action,Loglikel,opt] = pparrythmia(EKGR, varargin)
% function [EKGR,Action,Loglikel,opt] = pparrythmia(EKGR, varargin)
%
%
% Copyright (C) Luca Citi and Riccardo Barbieri, 2010-2011.
% All Rights Reserved. See LICENSE.TXT for license details.
% {lciti,barbieri}@neurostat.mit.edu
% http://users.neurostat.mit.edu/barbieri/pphrv

%% EXAMPLE
%data = load('Y2.pre');
%R = data(:,1);
%Re = R;
%Re = [Re(1:330); mean(Re(330:331)); Re(331:600); mean(Re(600:601)); Re(601:end)];
%Re([140 640]) = [];
%Re([230 440]) = Re([230 440])-.2;
%Re([280 680]) = Re([280 680])+.2;
%[Rc,Action,Loglikel,opt] = pparrythmia(Re);
%plot(R(2:end), diff(R))
%hold on
%plot(Re(2:end), diff(Re), 'r')
%plot(Rc(2:end), diff(Rc), 'g')

% Default options
opt.P = 5; % RR order
opt.hasTheta0 = 0; % wether or not the AR model has a theta0 constant to account for the average mu
opt.weight = 0.98; % weighting factor
opt.W = 60; % window length for local likelihood estimate (in seconds)
opt.fromFuture = 3;
opt.DL_add = [0 4];
opt.DL_remove = [3 8];
opt.DL_move = [2 7];
opt.DL_move2 = [8 28];
opt.DL_ignore = [6 14];
opt.DL_regr = [40 40];
opt.plainIgnore = 1;
opt.maximize_loglikel = @likel_invnorm_mex; % use loglikelihood of inverse gaussian
opt.debug = []; %%%%%%%%%%%%%%%

% THE FILE HAS BEEN CHANGED TO UNCOMMENT LINE 99 IN ORDER TO FIX THE GAPS on
% August 6 2013


% PROCESS OPTIONS
opt_names = fieldnames(opt);
J = 1;
while J <= length(varargin)
    key = varargin{J};
    keyN = find(strcmpi(key, opt_names));
    if isempty(keyN)
        warning('PPARRYTHMIA:UnknownOption', 'Do not know option ''%s''', key);
    else
        opt.(opt_names{keyN}) = varargin{J+1};
    end
    J = J + 2;
end

% proxies for options
W = opt.W;
P = opt.P;
maximize_loglikel = opt.maximize_loglikel;
thresholds = [opt.DL_add; opt.DL_remove; opt.DL_move; opt.DL_move2; opt.DL_ignore];
max_rr = 3; % if RR more than max_rr seconds, consider it a gap

opt.t0 = EKGR(1);
EKGR = EKGR(:) - opt.t0; % times are relative to EKGR(1) which is set to 0

lenR = length(EKGR);
j = find(EKGR(1:min(lenR,floor(5*W))) > W, 1) - 1; % index of last R event within W (check only first 5*W beats)
observ_ev = EKGR(1:j-1);
thetap = [];

steps = zeros(1, lenR);
Loglikel = NaN(7, lenR);
Action = zeros(1, lenR);

% try to initialize mask to avoid outliers
rr = diff(EKGR(1:j));
rr = abs(rr - median(rr));
Action(rr / median(rr) > 7) = 31;


tic();

while j <= length(EKGR)-2
    
    if mod(j, 1000) == 0
        fprintf(2, '%3ds: Processed %6d of %d\n', toc(), j, lenR);
    end
    time = EKGR(j);
    ii = find(EKGR(1:j) > time - W, 1);
    observ_ev = EKGR(ii:j);
    uk = observ_ev(P+2:end);
    rr = diff(observ_ev);
    if length(rr) < 2*P || max(rr) > max_rr
        if isempty(rr)
            Action(j) = 32; % gap
        end
        Action(find(rr > max_rr) - length(rr) + j) = 32; % gap
        j = j + 1;
        continue
    end
    wn = rr(P+1:end);
    xn = [];
    if opt.hasTheta0
        xn = ones(length(wn),1);
    end
    xn = [xn, toeplitz(rr(P:end-1), rr(P:-1:1))];
    eta = weights(time, uk, opt.weight);
    
    % mask for non-corrected outlying beats
    ll = sort(Loglikel(1,ii:j-1));
    ll = mean(max(ll, ll(round(.25 * length(ll)))));
    % mask bad regressors
    maskR = (Loglikel(1,ii:j) < ll - opt.DL_regr(2)) | (Action(ii:j) >= 16);
    mask = cumsum(maskR);
    mask = mask(P+2:end) == [0 mask(1:end-P-2)];
    % mask bad regressee (stricter)
    mask = mask & ~((Loglikel(1,ii+P+1:j) < ll - opt.DL_regr(1)) | (Action(ii+P+1:j) > 0));
    
    % check the likelihod of the new RR with the old model
    if ~isempty(thetap)
        if isnan(Loglikel(1,j)) && Action(j-1) == 0 && all(Action(j-[2 3]) ~= 4)
            %        if isnan(Loglikel(1,j)) && Action(j-1) < 16 && all(Action(j-(1:2)) ~= 4)
            xt = xn(end,:);
            if ~mask(end)
                last_P_ok = rr(find(~maskR(1:end-1), P, 'last')-1);
                xt(opt.hasTheta0+1:end) = max(min(last_P_ok), min(max(last_P_ok), xt(opt.hasTheta0+1:end)));
                %                badregr = maskR(end-1:-1:end-P);
                %                xt(opt.hasTheta0 + find(badregr)) = median(rr(find(~maskR(1:end-1), P, 'last')-1));
                %                badregr = Action(j-1:-1:j-P) >= 16;
                %                xt(opt.hasTheta0 + find(badregr)) = median(rr(find(~maskR(1:end-1), P, 'last')-1));
            end
            
            % fprintf('%.3f ', time)
            
            
            if ismember(j, opt.debug)
                [action, Rt1, loglikel] = arrythmia_invnorm_mex_debug(xt, EKGR(j-1:min(j+opt.fromFuture,length(EKGR))), opt.hasTheta0, thetap, k, thresholds);
            else
                [action, Rt1, loglikel] = arrythmia_invnorm_mex(xt, EKGR(j-1:min(j+opt.fromFuture,length(EKGR))), opt.hasTheta0, thetap, k, thresholds);
            end
            Loglikel(:,j) = loglikel([1 1 2 3 4 5 6]);
            if action > 0
                if opt.plainIgnore && (16 < action && action < 32)
                    action = 16;
                end
                switch action
                    case {1,17} % insert event
                        EKGR = [EKGR(1:j-1); Rt1; EKGR(j:end)];
                        lenR = lenR + 1;
                        Loglikel(:,end+1) = NaN;
                        Action(end+1) = 0;
                    case {2,18} % remove event
                        EKGR(j) = [];
                        lenR = lenR - 1;
                        Loglikel(:,end) = [];
                        Action(end) = [];
                    case {3,19} % move event
                        EKGR(j) = Rt1;
                    case {4,20} % move 2 events
                        EKGR(j:j+1) = Rt1;
                end
                Action(j) = action;
                %keyboard
                continue
            end
        else
            [tmp, tmp, tmp, tmp, loglikel] = maximize_loglikel(xn(end,:), wn(end), 1, thetap, k, 0*xn(end,:), 0, 0);
            Loglikel(1,j) = loglikel(1);
        end
    end
    
    if any(mask)
        [thetap, k, steps(j), loglikel] = maximize_loglikel(xn(mask,:), wn(mask), eta(mask));
        % fix AR models that might have become slightly unstable due to the estimation process
        % using an exponential decay (see Stoica and Moses, Signal Processing 26(1) 1992)
        %         poles = roots([1; -thetap]);
        %         mod_scale = min(.99/max(abs(poles)), 1);
        %         thetap = thetap .* (mod_scale .^ (1:length(thetap)))';
    end
    %    Thetap(:,j) = thetap;
    Kappa(j) = k;
    j = j + 1;
end

opt.Kappa = Kappa;

opt.steps = steps;
EKGR = EKGR + opt.t0;
%keyboard
