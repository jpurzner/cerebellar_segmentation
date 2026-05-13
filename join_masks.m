function junctions = join_masks(pia_deep_bwl, all_outline, pia_discard)

[outline_y, outline_x] = find(all_outline);  % x and y are column vectors.
outline_xy = [outline_y, outline_x];

pia_junction = zeros(size(pia_deep_bwl)); 
for n = 1:max(unique(pia_deep_bwl))
    if ~ismember(n,pia_discard)
        [~, ~, re, ce] = findendsjunctions(pia_deep_bwl == n);
        deep_ends = [re, ce];
        [D,I] = pdist2(outline_xy, deep_ends, 'euclidean', 'Smallest',1);
        [min_D, deep_end_I] = min(D);
        if min_D < 400 
            outline_closest = outline_xy(I(deep_end_I),:);
            deep_ends_closest = deep_ends(deep_end_I,:);
            conn_line = [outline_closest(2), outline_closest(1), deep_ends_closest(2), deep_ends_closest(1)];
            pia_junction = insertShape(double(pia_junction),'line',conn_line);
        end
    end
end 

junctions = (pia_junction(:,:,1) > 0);
return
end