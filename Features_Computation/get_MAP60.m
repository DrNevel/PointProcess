function [map60, missing, varargout]=get_MAP60(abp,Fs,varargin)

do_plot='off';
if nargin > 2
    for i= 1 : 1 : length(varargin)
        switch varargin{i}
            case 'plot'
                do_plot='on';
            case 'abp_start'
                abp_start = varargin{i+1}; % [samples]
        end
    end
end

if size(abp,1) > 1
    abp=abp';
end
if ~exist('abp_start','var'), abp_start = 1; end

abp = abp(abp_start:end);

T=Fs*60; % nr of samples in 1 minute
intrvl=floor(length(abp)/T); % nr of intervals (== nr of minutes)

t_map60=zeros(intrvl,1);
for i= 1 : 1 : intrvl
    t_map60(i)=i*T;
end

map60=find_map60(abp,length(abp),T,intrvl);

missing = find(isnan(map60));
if ~isempty(missing) && missing(1)-1 ~= 0 && missing(end) ~= intrvl
    map60(missing) = interp1([missing(1)-1,missing(end)+1]...
        ,[map60(missing(1)-1),map60(missing(end)+1)],missing);
end


if strcmp(do_plot,'on')
    figure, plot((1 : length(abp))./(60*125), abp, t_map60./(60*125),map60,'r');
end

if nargout==3
    varargout{1}=(t_map60(:)' + abp_start)/(125*60);
end

end
