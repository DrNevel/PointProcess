function [TOTPWR, ULF, VLF, LF, HF, LFtoHF,LFnu,HFnu, slope, out]=get_spect2(amp,ind,varargin)
% Every 5 minutes: 
% [TOTPWR, ULF, VLF, LF, HF, LFtoHF,LFnu,HFnu, slope, out]=get_spect2(amp,ind,'pyulear')
% Whole signal:
% [TOTPWR, ULF, VLF, LF, HF, LFtoHF,LFnu,HFnu, slope, out]=get_spect2(amp,ind,'plomb','all_sig')
%
% Last Update: 2021/02/12: AR order selection criteria
% Last Update: 2021/02/18: AR order selection criteria


warning('off')
opt=1;
all_sig=0;
t_diff = diff(ind);
out=[];
fixed_order=[];

if nargin > 2
    for i= 1 : 1 : length(varargin)
        switch varargin{i}
            case 'plomb'
                opt=0;
            case 'pyulear'
                opt=1;
            case 'all_sig'
                all_sig=1;
            case 'fixed_order'
                fixed_order=varargin{i+1};
        end
    end
end

if all_sig==0 % Compute every 5 minutes
    
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
                error('Plomb is not finalized, correct o change')
                [pxx{k,:}, f{k,:}]=plomb(sig,ind(start:i));
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
        orders = [7:2:13];
        cumulative = 0;
        start = 1;
        k = 1;
        for i = 1:length(t_diff)
            cumulative = cumulative + t_diff(i);
            if (cumulative >= window) || (cumulative>=window-60 && i==length(t_diff))
                cumulative = 0;
                sig = amp(start:i)-mean(amp(start:i));
                if size(sig(:),1)<order
                    continue;
                end
                if isempty(fixed_order)
                [AIC_v(k,:),FPE_v(k,:),PEWT(k,:)] = AR_assess(sig,orders);
                [~,opt_ord_1i] = min(AIC_v(k,:));
                opt_ord_1 = orders(opt_ord_1i) ;
                opt_ord_2 = min(orders(PEWT(k,:)==0));
                if isempty(opt_ord_2)
                    opt_ord(k)=opt_ord_1;
                elseif PEWT(k,opt_ord_1i)==0
                    opt_ord(k)=opt_ord_1;
                else
                    opt_ord(k)=opt_ord_2;
                end
                OPT.opt_ord_1=opt_ord_1;
                OPT.opt_ord_2=opt_ord_2;
                OPT.opt_ord(k)=opt_ord(k);
                [pxx(k,:), f(k,:)]=pyulear(sig,opt_ord(k),1024,1/mean(t_diff(start:i),'omitnan'));
                else
                    [pxx(k,:), f(k,:)]=pyulear(sig,fixed_order,1024,1/mean(t_diff(start:i),'omitnan'));
                    OPT.fixed_order=fixed_order;
                end
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
            out.fitting.PEWT = PEWT;
            out.fitting.opt_ord = opt_ord;
            
            Pxx = mean(pxx,1);
            F = f(end,:);
        else
            %         disp('o');
        end
    end
    
elseif all_sig==1 % Use the whole signal
    
    if opt==0
        sig = amp-mean(amp);
        [Pxx, F]=plomb(sig,ind);
        out.Pxx=Pxx;
        out.F=F;
    else
        order = 9;
        sig = amp-mean(amp);
        if isempty(fixed_order)
        orders = [7:15];
        [AIC_v,FPE_v,PEWT] = AR_assess(sig,orders);
        [~,opt_ord_1i] = min(AIC_v);
        opt_ord_1 = orders(opt_ord_1i) ;
        opt_ord_2 = min(orders(PEWT==0));
        if isempty(opt_ord_2)
            opt_ord=opt_ord_1;
        elseif PEWT(opt_ord_1i)==0
            opt_ord=opt_ord_1;
        else
            opt_ord=opt_ord_2;
        end
        OPT.opt_ord_1=opt_ord_1;
        OPT.opt_ord_2=opt_ord_2;
        OPT.opt_ord=opt_ord;
        [Pxx, F]=pyulear(sig,opt_ord,1024,1/mean(t_diff));
        else
        [Pxx, F]=pyulear(sig,fixed_order,1024,1/mean(t_diff));
        OPT.fixed_order=fixed_order;
        end
        out.Pxx=Pxx;
        out.F=F;
    end
end

if exist('OPT','var')
    out.OPT=OPT;
end
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5624990/
if exist('Pxx','var')
    f=diff(F(F <= 0.4));
    if size(f,1)>=size(f,2)
        f=[f(1);f];
    else
        f=[f(1),f];
    end
    TOTPWR=sum(Pxx(F <= 0.4).*f);
    f=diff(F(F <= 0.003));
    if size(f,1)>=size(f,2)
        f=[f(1);f];
    else
        f=[f(1),f];
    end
    ULF=sum(Pxx(F <= 0.003).*f);
    f=diff(F(F > 0.003 & F <= 0.04));
    if size(f,1)>=size(f,2)
        f=[f(1);f];
    else
        f=[f(1),f];
    end
    VLF=sum(Pxx(F > 0.003 & F <= 0.04).*f);
    f=diff(F(F > 0.04 & F <= 0.15));
    if size(f,1)>=size(f,2)
        f=[f(1);f];
    else
        f=[f(1),f];
    end
    LF=sum(Pxx(F > 0.04 & F <= 0.15).*f);
    f=diff(F(F > 0.15 & F <= 0.4));
    if size(f,1)>size(f,2)
        f=[f(1);f];
    else
        f=[f(1),f];
    end
    HF=sum(Pxx(F > 0.15 & F <= 0.4).*f);
    LFtoHF=LF/HF;
    % LFnu = LF/(LF+HF);
    % HFnu = HF/(LF+HF);
    LFnu = LF/(TOTPWR-VLF);
    HFnu = HF/(TOTPWR-VLF);
    
    % Use all Info available with a LOMB Periodogram
    [SP, FR]=plomb(amp,ind);
    y = SP(FR <= 0.04);
    f_2 = FR(FR <= 0.04);
    slope = polyfit(log(f_2(2:end)),log(y(2:end)),1);
    slope = slope(1);
else
    [TOTPWR, ULF, VLF, LF, HF, LFtoHF,LFnu,HFnu, slope, out] = deal([]);
end
% To verify:
% Y = polyval(slope,log(f_2));
% plot(log(f_2),log(y),log(f_2),Y)

warning('on')
end







