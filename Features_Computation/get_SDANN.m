function [ SDANN ] = get_SDANN( hrv )
% Standard deviation of the averages of NN intervals in all 5-minute 
% segments of the entire reconding.
% hrv in seconds.

i = 1;
AVRGS = [];

for f = 2:length(hrv)
    time = sum(hrv(i:f));
    if time >= 5*60
        AVRGS = cat(1,AVRGS,mean(hrv(i:f)));
        i = f+1;
    end
    if f+1 == length(hrv)
        AVRGS = cat(1,AVRGS,mean(hrv(i:f+1)));
        break;
    end
end

SDANN = std(AVRGS)*1000;

end

