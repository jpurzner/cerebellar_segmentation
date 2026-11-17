function [ im_new ] = keep_largest( im )
%keeps the largest connected component
%   Detailed explanation goes here

im_bwl = bwlabel(im);
im_bwl_counts = histcounts(im_bwl, max(unique(im_bwl)));
[~, im_max_lab] = max(im_bwl_counts(2:end));
im_new = im_bwl == im_max_lab;
return

end

