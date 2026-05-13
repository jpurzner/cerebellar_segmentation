function [ bw_ov ] = im_get_overlap(bw1, bw2, ov_tresh, gap )
%im_get_intersected for a mask get the overlaping segments 

% ov_thres is the overlap threshold that is used to accept part of bw2
% gap is the size between grid used to chop up bw2 

% bw_ov is the sum of overalps 

% make a cut around the image using the perimeter of bw1
% this will make more segments and split unwanted regions 
%bw1_perim = bwperim(bw1);
%bw2 = bw2 - bw1_perim;
%bw2 = bw2 == 1;


% weird solution to the problem of having spurious connections 
% make a grid of 0s that move over the image, this breaks up bw2 to prevent
% over inclusive overlaps
bw_ov_t = zeros(size(bw1));

for start = 1:(gap/10):gap 
    bw2_t = bw2;
    % cut bw2 into pieces using a regular grid 
    for n = start:gap:size(bw2_t,1)
       bw2_t(n,:) = 0; 
    end    
    for n = start:gap:size(bw2_t,2)
       bw2_t(:,n) = 0; 
    end    

    bw2_bwl = bwlabel(bw2_t);
    ov_idx = unique(bw1 .* bw2_bwl);
    ov_idx = sort(ov_idx);
    ov_idx = ov_idx(2:end);

    ov_idx_quant = histc(reshape(bw2_bwl .* bw1,[],1), ov_idx) ./ histc(reshape(bw2_bwl ,[],1), ov_idx);

    % keep high overlaping CSF space 
    ov_sel = ov_idx(ov_idx_quant > ov_tresh);

    bw_ov_t = bw_ov_t + ismember(bw2_bwl, ov_sel );
    %bw_ov = ismember(bw2_bwl , ov_idx(2:end));
end 

bw_ov = bw_ov_t;

return

end

