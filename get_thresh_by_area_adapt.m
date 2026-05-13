function [ sens ] = get_thresh_by_area_adapt( im, keep_mask, area_frac )
%thresh_by_area thresholds the image by the area of a mask that is marked  

%high_i = min(find(cumsum(pixelCounts)/sum(sum(keep_mask)) > (1- area_frac)));
%high_in = grayLevels(high_i+1);

% set the values outside the mask to NaN 
im = double(im);
im_mask = im;
im_mask(~keep_mask) = nan;
im_mask =   im_mask .*  keep_mask; 

sensitivities = 0.02:0.02:0.9;
area = zeros(size(0.1:0.05:0.9));
count = 1;
for s = sensitivities
    %disp(s)
    T = adaptthresh(im , s,'NeighborhoodSize' ,  2*floor(size(im)/25)+1);
    bin = imbinarize(im ,T);
    area(count) = sum(sum(bin .* keep_mask))/sum(sum(keep_mask));
    count = count + 1;
end

[ d, ix ] = min( abs( area-(area_frac) ) );
sens = sensitivities(ix);
return
%high_i = min(find(cumsum(pixelCounts)/sum(sum(mask)) > (1- area_frac)));
%high_in = grayLevels(high_i+1);
%im_bin = im >= high_in;
return
end

