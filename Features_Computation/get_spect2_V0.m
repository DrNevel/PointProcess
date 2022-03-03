function [TOTPWR, ULF, VLF, LF, HF, LFtoHF,LFnu,HFnu, slope, out]=get_spect2(amp,ind,varargin)

opt=1;
t_diff = diff(ind);

if nargin > 2
    for i= 1 : 1 : length(varargin)
        switch varargin{i}
            case 'plomb'
                opt=0;
            case 'pyulear'
                opt=1;
        end
    end
end

window = 5*60; % 5 minutes

if opt==0
    while sum(diff(ind)<=0)~=0
        amp(diff(ind)<=0)=[];
        ind(diff(ind) <=0)=[];
    end
    cumulative = 0;
    start = 1;
    k = 1;
    for i = 1:length(t_diff)
        cumulative = cumulative + t_diff(i);
        if cumulative >= window
            cumulative = 0;
            sig = amp(start:i)-mean(amp(start:i));
            [pxx(k,:), f(k,:)]=plomb(sig,ind(start:i));
            start = i+1;
            k = k+1;
        elseif i == length(amp)
            break;
        end
    end
    Pxx = mean(pxx,1);
    F = f(end,:);
else
    order = 9;
    orders = [5:2:15];
    cumulative = 0;
    start = 1;
    k = 1;
    for i = 1:length(t_diff)
        cumulative = cumulative + t_diff(i);
        if cumulative >= window
            cumulative = 0;
            sig = amp(start:i)-mean(amp(start:i));
            if size(sig(:),1)<order
                continue;
            end
            [AIC_v(k,:),FPE_v(k,:),PEWT(k,:)] = AR_assess(sig,orders);
            [~,opt_ord] = min(AIC_v(k,:));
            [pxx(k,:), f(k,:)]=pyulear(sig,orders(opt_ord),1024,1/mean(t_diff(start:i),'omitnan'));
            start = i+1;
            k = k+1;
        elseif i == length(amp)
            break;
        end
    end
    % Time-varying powers
    if exist('pxx','var')
        out.TOTPWR = sum(pxx(:,f(end,:) <= 0.4),2);
        out.ULF = sum(pxx(:,f(end,:) <= 0.003),2);
        out.VLF = sum(pxx(:,f(end,:) > 0.003 & f(end,:) <= 0.04),2);
        out.LF = sum(pxx(:,f(end,:) > 0.04 & f(end,:) <= 0.15),2);
        out.HF = sum(pxx(:,f(end,:) > 0.15 & f(end,:) <= 0.4),2);
        out.LFtoHF = out.LF./out.HF;
        out.fitting.orders = orders;
        out.fitting.AIC = AIC_v;
        out.fitting.FPE = FPE_v;
        
        Pxx = mean(pxx,1);
        F = f(end,:);
    else
%         disp('o');
    end
end

% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5624990/
if exist('Pxx','var')
    
    TOTPWR=sum(Pxx(F <= 0.4));
    ULF=sum(Pxx(F <= 0.003));
    VLF=sum(Pxx(F > 0.003 & F <= 0.04));
    LF=sum(Pxx(F > 0.04 & F <= 0.15));
    HF=sum(Pxx(F > 0.15 & F <= 0.4));
    LFtoHF=LF/HF;
    % LFnu = LF/(LF+HF);
    % HFnu = HF/(LF+HF);
    LFnu = LF/(TOTPWR-VLF);
    HFnu = HF/(TOTPWR-VLF);
    
    y = Pxx(F <= 0.04);
    f_2 = F(F <= 0.04);
    slope = polyfit(log(f_2(2:end)),log(y(2:end)),1);
    slope = slope(1);
else
    [TOTPWR, ULF, VLF, LF, HF, LFtoHF,LFnu,HFnu, slope, out] = deal([]);
end
% To verify:
% Y = polyval(slope,log(f_2));
% plot(log(f_2),log(y),log(f_2),Y)
end







