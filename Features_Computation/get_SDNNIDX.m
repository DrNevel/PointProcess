function [ SDNNIDX ] = get_SDNNIDX( hrv )
% Mean of the standard deviations of NN intervals in all 5-minute segments 
% of the entire recording
% hrv in seconds.

i = 1;
SD = [];

for f = 2:length(hrv)
    time = sum(hrv(i:f)); % Times are in seconds
    if time >= 5*60
        SD = cat(1,SD,std(hrv(i:f)));
        i = f+1;
    end
    if f+1 == length(hrv)
        SD = cat(1,SD,std(hrv(i:f+1)));
        break;
    end
end

SDNNIDX = mean(SD)*1000;

end

