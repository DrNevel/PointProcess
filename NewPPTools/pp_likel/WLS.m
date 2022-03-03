function [thetap,var] = WLS(xn, wn, p)
% weighted least-squares regression
%
% Copyright (C) Maximiliano Mollura and Riccardo Barbieri, 2019-2020.
% All Rights Reserved. See LICENSE.TXT for license details.
% maximiliano.mollura@polimi.it
% riccardo.barbieri@polimi.it

if ~exist('p','var')
    p = ones(1,size(wn,1));
else
    assert(length(p)==size(wn,1));
end

P = diag(p);

thetap = (xn'*P*xn)\(xn'*P*wn);
var = sum(p(:).*(wn - xn*thetap).^2)/sum(p);