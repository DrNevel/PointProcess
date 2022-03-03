function results=BIVARIATE_AR(RR,RESP,ord_rr,ord_resp)
% Ax=b -> A*OBS=RR_T_1
% RR(i)=RR(i+1)-RR(i)
% RESP(i)

if nargin<1
    RR=1:999;
    RESP=[1:1000].*10;
    ord_rr=9;
    ord_resp=9;
end
RR_T=toeplitz(RR(ord_rr:end-1),RR(ord_rr:-1:1));
RR_T_1=RR(ord_rr+1:end)';
RESP_T=toeplitz(RESP(ord_resp:end-2),RESP(ord_resp:-1:1));
RESP_T_1=RESP(ord_rr+1:end-1)';

OBS=[ones(size(RR_T,1),1),RR_T,RESP_T];

coeffs_RESPtoRR=(OBS'*OBS)\OBS'*RR_T_1;

coeffs_RRtoRESP=(OBS'*OBS)\OBS'*RESP_T_1;

results.TOEPLITZ=OBS;
results.out_rr=RR_T_1;
results.out_rest=RESP_T_1;
results.coeffs_RESPtoRR=coeffs_RESPtoRR;
results.coeffs_RRtoRESP=coeffs_RRtoRESP;

end