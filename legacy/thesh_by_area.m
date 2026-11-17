function [ im_bin ] = thesh_by_area( im, mask, area_frac )
%thresh_by_area thresholds the image by the area of a mask that is marked  

[pixelCounts , grayLevels] = imhist(im,200);
[max_pixelCounts, max_i] =  max(pixelCounts);
% get high point of histogram (noise) and add the typical width 
% clip the low pixel counts 
high_i = min(find(cumsum(pixelCounts(low_i:end))/sum(sum(mask))) > area_frac));
high_i = low_i + high_i;
high_in = grayLevels(high_i);
im_bin = im_bin >= high_in;
return
end

