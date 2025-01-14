function [peak_final, location_final] = Peak_find(De_cond_hilbert, Threshold, slide_window)
% Function£ºFind the peak and its corresponding position based on the envelope of the Hilbert transform.
% Input£º
% De_cond_hilbert - MAE signal after Hilbert transform
% Threshold -
% slide_window - By setting a sliding window, we believe that there will not be two conductivity changes within a narrow range.
%                In 2D configuration, the segmentation is adopted.
% Ouput£º
% peak_final - Found peak
% location_final - Found location
% Author£ºTong Sun
% Date£º2022-9-23

[peaks, location] = findpeaks(De_cond_hilbert);     % Find the peak from the data after Hilbert transform, where all peaks and positions will be found
location_thres = location(peaks > Threshold);       % Those below the threshold will be excluded, and the position will be returned here
peaks_thres = peaks(peaks > Threshold);             % Those below the threshold will be excluded, and the peak value will be returned here
peak_final = [];                                    % Peak storage
location_final = [];                                % Location storage
i = 1;                                              % Loop variable
location_temp = location_thres;                     % Location temporary variable

while(i <= size(location_thres,1))
    index = (location_thres(i,1)+slide_window) >= location_temp;
    temp = location_temp(index);
    peak_final = [peak_final;max(De_cond_hilbert(temp))];   
    location_final = [location_final; location_thres(peaks_thres==max(De_cond_hilbert(temp)))]; 
    location_temp(1:size(temp,1)) = [];
    i = i + size(temp,1);
end

end