function [ high_in ] = get_thresh_by_area( im, keep_mask, area_frac )
%thresh_by_area thresholds the image by the area of a mask that is marked  

% set the values outside the mask to NaN 
im = double(im);
im(~keep_mask) = nan;
im =   im .*  keep_mask; 
[pixelCounts , grayLevels] = histcounts(im,1000);
% get high point of histogram (noise) and add the typical width 
% clip the low pixel counts 

high_i = min(find(cumsum(pixelCounts)/sum(sum(keep_mask)) > (1- area_frac)));
high_in = grayLevels(high_i+1);

return
end

