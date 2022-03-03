function [NN50, pNN50]=get_NN50(hrv)
% NN50: Number of successive RR interval pairs that differ more than 50 ms
% pNN50: NN50 divided by the total number of RR intervals [%]
% See: http://bsamig.uef.fi/kubios/kubios_hrv_users_guide.pdf


diff_hrv=diff(hrv);

spot_nn50=zeros(1,length(diff_hrv));
spot_nn50(abs(diff_hrv) > 0.05)=1;

NN50=sum(spot_nn50);

pNN50=(NN50*100)/(length(hrv)-1);

end