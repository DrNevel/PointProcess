function [ SDNN ] = get_SDNN( hrv )
% hrv in seconds.

SDNN = std(hrv)*1000;

end

