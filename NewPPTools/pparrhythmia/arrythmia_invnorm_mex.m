function varargout = arrythmia_invnorm_mex(varargin)
% [action, Rt1, loglikel] = arrythmia_invnorm_mex(xt, R(j-1:j+1), hasTheta0, thetap, k, thresholds);
%
%
% Copyright (C) Luca Citi and Riccardo Barbieri, 2010-2011.
% All Rights Reserved. See LICENSE.TXT for license details.
% {lciti,barbieri}@neurostat.mit.edu
% http://users.neurostat.mit.edu/barbieri/pphrv

% if we get here it means the mex file does not exist,
% we need to build it

if exist('arrythmia_invnorm_mex', 'file') ~= 3
    curpath = pwd();
    fpath = fileparts(which('arrythmia_invnorm_mex'));
    cd(fpath);
    mex('arrythmia_invnorm_mex.c');
    cd(curpath);
    rehash;
end

if exist('arrythmia_invnorm_mex', 'file') == 3
    varargout = cell(max(nargout, 1), 1);
    [varargout{:}] = arrythmia_invnorm_mex(varargin{:});
end
