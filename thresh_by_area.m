function [ im_bin ] = thresh_by_area( im, mask, area_frac )
%thresh_by_area thresholds the image by the area of a mask that is marked  

% set the values outside the mask to NaN 
im = double(im);
im(~(mask)) = nan;
[pixelCounts , grayLevels] = histcounts(im,500);
% get high point of histogram (noise) and add the typical width 
% clip the low pixel counts 

high_i = min(find(cumsum(pixelCounts)/sum(sum(mask)) > (1- area_frac)));
high_in = grayLevels(high_i+1);
im_bin = im >= high_in;
return
end

