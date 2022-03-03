function SDSD=get_SDSD(hrv)
% Standard Deviation of the differences of the NN intervals

d_hrv=diff(hrv);
SDSD=std(d_hrv)*1000;
end