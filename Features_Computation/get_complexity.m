function [D,Alpha1,H,ApEn,SampEn,CorrDim,LyapExp,XEntropy,MSEntr] = get_complexity(amp, dimension, lag, radius,amp_sig2)
fs = 125;
%% Aggregate Series
% try
%     Agg4 = aggregate_series(amp,4); % 15 mins
% catch
%     Agg4 = NaN(1,4);
% end
% try
%     Agg5 = aggregate_series(amp,5); % 12 mins
% catch
%     Agg5 = NaN(1,5);
% end
% try
%     Agg6 = aggregate_series(amp,6); % 10 mins
% catch
%     Agg6 = NaN(1,6);
% end
%% Detrended Fluctuation Analysis
try
    [D,Alpha1]=DFA_main(amp);
catch
    D = NaN(1,1);
    Alpha1 = NaN(1,1);
end
%% Hurst Exponent
try
    H = higuchi(amp,0);
catch
    H = NaN(1,1);
end
%% Entropy Measures
try
    if size(amp(:),1) > 50000
        amp_2 = amp(1:4:end);
    else
        amp_2 = amp;
    end
    ApEn = approximateEntropy(amp_2, 'Dimension', dimension, 'Lag', lag, 'Radius', radius); % MATLAB FCN
catch
    ApEn = NaN(1,1);
end
try
    if size(amp(:),1) > 50000
        amp_2 = amp(1:4:end);
    else
        amp_2 = amp;
    end
    SampEn = sampen(amp_2,dimension,radius,'euclidean'); % https://it.mathworks.com/matlabcentral/fileexchange/35784-sample-entropy
catch
    SampEn = NaN(1,1);
end
%% Correlation Dimension
try
    if size(amp(:),1) > 50000
        amp_2 = amp(1:4:end);
    else
        amp_2 = amp;
    end
    CorrDim = correlationDimension(amp_2, 'Dimension', dimension, 'Lag', lag, 'MaxRadius', radius);
catch
    CorrDim = NaN(1,1);
end
%% Lyapunov Exponent
% LyapExp = NaN(1,1);
try
    if size(amp(:),1) > 50000
        amp_2 = amp(1:4:end);
    else
        amp_2 = amp;
    end
    LyapExp = lyapunovExponent(amp_2, 1/mean(amp,'omitnan'), 'Dimension', dimension, 'Lag', lag);
catch
    LyapExp = NaN(1,1);
end
%% Cross-Sample Entropy
XEntropy = NaN(1,1);
if exist('amp_sig2','var')
    try
        if size(amp(:),1) > 50000
            amp_2 = amp(1:4:end);
            amp_sig2 = amp_sig2(1:4:end);
        else
            amp_2 = amp;
            amp_sig2=amp_sig2;
        end
        XEntropy = cross_sampen(amp_2,amp_sig2,dimension,radius,1);
    catch
        XEntropy = NaN(1,1);
    end
end
%% Multiscale Entropy - It Needs very long samples otherwise no differences can be found
try
    MSEntr = msentropy(amp_2,[],[],[],[],[],[],[],10);
    MSEntr = MSEntr(:)';
catch
    MSEntr = NaN(1,10);
end



