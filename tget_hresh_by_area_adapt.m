function [ im_bin ] = thresh_by_area_adapt( im, mask, area_frac )
%thresh_by_area thresholds the image by the area of a mask that is marked  

% set the values outside the mask to NaN 
im(~(mask)) = nan;

area = zeros(size(0.1:0.05:0.9));
count = 1;
for s = 0.1:0.05:0.9
    T = adaptthresh(im_bin , s,'NeighborhoodSize' , 2*floor(size(im)/35)+1);
    bin = imbinarize(im_bin ,T);
    area(count) = sum(sum(bin))/sum(sum(mask));
end

im_bin = area;

%high_i = min(find(cumsum(pixelCounts)/sum(sum(mask)) > (1- area_frac)));
%high_in = grayLevels(high_i+1);
%im_bin = im >= high_in;
return
end

