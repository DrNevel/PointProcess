function [NN20, pNN20]=get_NN20(hrv)
% NN20: Number of successive RR interval pairs that differ more than 20 ms
% pNN20: NN20 divided by the total number of RR intervals [%]
% See: http://bsamig.uef.fi/kubios/kubios_hrv_users_guide.pdf

diff_hrv=diff(hrv);

spot_nn20=zeros(1,length(diff_hrv));
spot_nn20(abs(diff_hrv) > 0.02)=1;

NN20=sum(spot_nn20);

pNN20=(NN20*100)/(length(hrv)-1);

end