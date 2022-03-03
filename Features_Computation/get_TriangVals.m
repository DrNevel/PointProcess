function [TRI,TINN] = get_TriangVals(hrv)
% Baseline width of the RR interval histogram
% See: https://github.com/MarcusVollmer/HRV/blob/master/HRV.m

% hrv = TaskForce.EKGR{1,1}  

[TRI,TINN]  = HRV.triangular_val(hrv,0,1/128,0);
TINN = TINN*1000;

end