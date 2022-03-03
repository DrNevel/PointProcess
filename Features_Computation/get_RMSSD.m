function RMSSD=get_RMSSD(hrv)
% Square root of the mean squared differences between successive RR [ms]
% intervals. See: http://bsamig.uef.fi/kubios/kubios_hrv_users_guide.pdf

d2=hrv(2:end);
d1=hrv(1:end-1);
d2_1=d2-d1;
RMSSD = rms(d2_1)*1000;

end