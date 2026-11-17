function layer = cerebellum_superpixels2(file)  

a = imread(file,1);
b = imread(file, 2);
c = imread(file, 3);

ha = homomorphic(a, 2, 0.25,2,0,5);
hb = homomorphic(b, 2, 0.25,2,0,5);
%hc = homomorphic(c, 2, 0.25,2,0,5);

ha_g =  imgaussfilt(ha,5);
hb_g =  imgaussfilt(hb,5);
%hc_g =  imgaussfilt(hc,5);

%col_g(:,:,1) = ha_g;
%col_g(:,:,2) = hb_g;
%col_g(:,:,3) = hc_g;

%[L,N] = superpixels(col_g,20000);


a_n = a;
a_n(a_n< 8000) =NaN;

b_n = b;
b_n(b_n< 8000) = NaN;


c_n = c;
c_n(c_n< 8000) = NaN;


a_norm = homomorphic(a_n, 2, 0.25,2,0,5);
%a_norm = (double(a_n)./double(b_n)) + double(a_n) ;
%a_norm(a_norm == Inf) = NaN;
%a_norm(a_norm == -Inf) = NaN;

c_norm = (double(c_n)./double(b_n)) + (2*double(c_n)) ;
c_norm(c_norm == Inf) = NaN;
c_norm(c_norm == -Inf) = NaN;


% totally flatten the images prior to division 
%hb_n = homomorphic(b_n, 2, 0.25,2,0,5);
%ha = homomorphic(a, 2, 0.25,2,0,5);
b_icor = illumination_cor(b);
% this part also removes grime or bleach marks 
a_nf = adapthisteq(mat2gray(a),'clipLimit',0.01,'Distribution','rayleigh', 'NumTiles', [50 50], 'Range', 'original');
b_nf = adapthisteq(mat2gray(b_icor),'clipLimit',0.01,'Distribution','rayleigh', 'NumTiles', [50 50], 'Range', 'original');
a_nfg =  imgaussfilt(a_nf,5);
b_nfg =  imgaussfilt(b_nf,5);

% remove low intensity parts of image 
a_nfg = mat2gray(a_nfg, [0.4 1]); 
b_nfg = mat2gray(b_nfg, [0.4 1]);

b_norm = (double(b_nfg) - double(a_nfg))   ;


%b_norm = (double(hb_n)./double(ha))   ;
%b_norm(b_norm == Inf) = NaN;
%b_norm(b_norm == -Inf) = NaN;



%b_norm(isnan(b_norm))=0;
%b_norm = b_norm.^2;
%b_norm = ordfilt2(b_norm,16,ones(4,4));
%b_norm =  imgaussfilt(b_norm,4);
%b_norm = illumination_cor(b_norm);
%b_norm = mat2gray(b_norm);


c_norm(isnan(c_norm))=0;
c_norm = c_norm.^2;
c_norm = ordfilt2(c_norm,16,ones(4,4));
c_norm = imgaussfilt(c_norm,5);
c_norm = mat2gray(c_norm);


bb = mat2gray(b);
all_cerebellum = bb > 0.1;
%all_cerebellum = imclearborder(all_cerebellum, 4);
all_cerebellum = imfill(all_cerebellum, 'holes');


% seperate the deep white matter from the molecular layer 
ha_n = homomorphic(a_n, 2, 0.25,2,0,5);
% this seem trivial but you can't have too many zeros in the denominator 
% or the ratio image looks weird
hc = homomorphic(c, 2, 0.25,2,0,5);
dwl_norm = log2(ha_n./hc);
dwl_norm = anisodiff(dwl_norm, 15, 50, 0.1, 1) ;



dwl_norm = mat2gray(dwl_norm);

% get the molecular layer
hb = homomorphic(b, 2, 0.25,2);
ha = homomorphic(a, 2, 0.25,2);
ml_n = hb .* ha;



%ml = (b/1000).^2;
%ml_n = anisodiff(ml_n, 15, 75, 0.1, 2);
ml_n =  imgaussfilt(ml_n,6);
ml_n =  ml_n.^2;
ml_n =  mat2gray(ml_n);
ml_n = 1 - ml_n;
ml_n = ml_n .* all_cerebellum;




% crudely segregate iEGL
%a_norm = a_norm .* (1 - double(c_bin));
a_norm = anisodiff(a_norm, 15, 75, 0.1, 2);

a_norm(isnan(a_norm))=0;


%a_norm = a_norm.^2;
%a_norm = ordfilt2(a_norm,16,ones(4,4));
%a_norm = imgaussfilt(a_norm,5);
a_norm = mat2gray(a_norm);


col_g(:,:,1) = mat2gray(a_norm,[0.5 1]) .* all_cerebellum ;
col_g(:,:,2) = mat2gray(b_norm, [0.02 0.25]) .* all_cerebellum;

col_g(:,:,2) = col_g(:,:,2) .* all_cerebellum;
col_g(:,:,3) = mat2gray(c_norm, [0.1 0.3]) .* all_cerebellum;

%layer = col_g; 
%return 

[L,N] = superpixels(col_g,20000);
sp_vals = regionprops(L,'PixelIdxList','PixelList', 'centroid');


image_cell = {col_g(:,:,1), col_g(:,:,2), col_g(:,:,3), mat2gray(ml_n, [0.4 0.75]), mat2gray(dwl_norm, [0.8 0.95]), double(all_cerebellum)}; 


sum_im_cell = repmat({zeros(1,N)},1,size(image_cell,2));
median_im_cell = repmat({zeros(1,N)},1,size(image_cell,2));
mean_im_cell = repmat({zeros(1,N)},1,size(image_cell,2));

% need to vectorize this nested loop 
for n=1:length(sp_vals)  %for each cell in the image
    for i=1:size(image_cell, 2)
        sum_im_cell{i}(n) = sum(image_cell{i}(sp_vals(n).PixelIdxList));
        median_im_cell{i}(n) = median(image_cell{i}(sp_vals(n).PixelIdxList));
        mean_im_cell{i}(n) = mean(image_cell{i}(sp_vals(n).PixelIdxList));
    end 
end


%sum_struct = cell2struct(sum_im_cell, strcat(names, '_sum'),2);
%median_struct = cell2struct(median_im_cell, strcat(names, '_median'),2);
%sp_vals = catstruct(sp_vals, median_struct);
%sp_vals = catstruct(sp_vals, sum_struct);


median_im = cell2mat(median_im_cell');
mean_im = cell2mat(mean_im_cell');
mean_im = mean_im';


%M = containers.Map(1:length(k),k, 'UniformValues',false, 'KeyType', 'uint16', 'ValueType', 'uint16');
%M = containers.Map('KeyType', 'uint32', 'ValueType', 'uint16');
%M = containers.Map(uint32(1:length(k)),uint16(k), 'UniformValues',true);

% start with heuristics to get the best matching superpixels then fill in
% remainder by proximity

%outer EGL
b_sp = find(median_im(2,:) >= 0.2 &  median_im(1,:) <  0.5);

%inner EGL
a_sp = find(median_im(1,:) >= 0.5 &  median_im(3,:) == 0);

%IGL 
c_sp = find(median_im(3,:) ~= 0 &  median_im(2,:) <  0.4);

%ML 
ml_sp =  find(median_im(4,:) > 0.7 & median_im(1,:) < 0.5 &  median_im(1,:) < 0.5);

%Deep white layer
dwl_sp =  find(median_im(5,:) > 0.2 & median_im(1,:) < 0.4 );

% cerebellum mask in superpixels 
cbl_sp = find(median_im(6,:) == 1);



%re-lable the superpixels 
sp_all = zeros(size(col_g,1),size(col_g,2));
sp_lab_c = cellstr(repmat('none',length(sp_vals),1));
sp_all(ismember(L, cbl_sp)) = 1;
sp_lab_c(cbl_sp) = cellstr('unmarked');
sp_all(ismember(L, ml_sp)) = 5;
sp_lab_c(ml_sp) = cellstr('ML');
sp_all(ismember(L, dwl_sp)) = 6;
sp_lab_c(dwl_sp) = cellstr('DWL');
sp_all(ismember(L, c_sp)) = 4;
sp_lab_c(c_sp) = cellstr('IGL');
sp_all(ismember(L, b_sp)) = 3;
sp_lab_c(b_sp) = cellstr('oEGL');
sp_all(ismember(L, a_sp)) = 2;
sp_lab_c(a_sp) = cellstr('iEGL');

% create a mapping from bwlabel to structure 
keySet =   [0,1,2,3,4,5,6];
valueSet = {'none','unmarked', 'iEGL','oEGL','IGL','ML','DWL'};
bwlab2layer_map = containers.Map(keySet,valueSet);
layer2bwlab_map = containers.Map(valueSet,keySet);

% skeletonize the outer egl to get a the outer sufcace 

all_outline = bwmorph(all_cerebellum,'remove');
se = strel('disk',20);
all_outline_mask= imclose(all_outline, se);
all_outline_mask = bwmorph(all_outline_mask,'thicken', 50);
all_outline_mask = bwmorph(all_outline_mask,'bridge', 50);
se = strel('disk',5);
all_outline_mask= imclose(all_outline_mask, se);



pia = imbinarize(b_norm, 'adaptive', 'Sensitivity', 0.6);
pia = bwareaopen(pia,200);
%se = strel('disk',30);
%pia= imclose(pia, se);

folia_skel = bwmorph(pia,'skel',Inf);
folia_skel = bwareaopen(folia_skel,200);
folia_skel = folia_skel .* (1-all_outline_mask);

pia_outline = all_outline + folia_skel;
pia_outline = bwmorph(pia_outline,'skel',Inf);
pia_outline = bwareaopen(pia_outline,500);
    
[edgelist edgeim, etype] = edgelink(pia_outline == 1);
tt = ismember(edgeim, find(etype == 1));



sp_all(pia_outline == 1) = 7;

imagesc(sp_all)

% assign the unmapped superpixels
sp_centroids =reshape([sp_vals.Centroid], 2, length(sp_vals));
sp_centroids = sp_centroids';

[pia_y,pia_x] = find(pia_outline);
pia_xy = [pia_x, pia_y];

% get minimum distance to pia for every centroid
pia_dist = pdist2( pia_xy, sp_centroids,'euclidean' ,'Smallest',1);
pia_dist = pia_dist';

% summarise the already assigned superpixels to get a center for the IF for
% each layer, centers are a vector of 5 values 
% AND pdist 
sp_t = array2table(median_im');

sp_t.pia_dist = pia_dist;
sp_t.label = sp_lab_c;
layer_centers = grpstats(sp_t,'label');


% NOTE for the distances the index starts with the iEGL so at the third
% entr
layer_centers_if = layer_centers{valueSet(3:7),3:7};
layer_centers_pia_dist = layer_centers{valueSet(3:7),'mean_pia_dist'};

 
pia_layer_dist = pdist2( layer_centers_pia_dist, pia_dist);
if_dist = pdist2(layer_centers_if, median_im(1:5,:)');

[min_pia_dist,min_pia_i] = min(pia_layer_dist',[],2);
[min_if_dist,min_if_i] = min(if_dist',[],2);


% reassign unmarked 
% by IF signal
change_idx = find(string(sp_lab_c) == 'unmarked' &  min_if_dist <=0.4);
change_label = horzcat(change_idx ,min_if_i(change_idx));
% by proximity to pia 
change_idx = find(string(sp_lab_c) == 'unmarked' &  min_pia_dist <=40);
change_label = vertcat(change_label, horzcat(change_idx, min_pia_i(change_idx)));
% ensure the pia mapping doesn't totally wreck the IF signal
if_new_d = arrayfun(@(bw_i, label_i) if_dist(label_i, bw_i), change_label(:,1), change_label(:,2));
change_label = change_label(if_new_d <= 0.5,:);


% specifically reassign the DWL classification that are within the EGL or
% ML (2,3,5) in sp_all, these sp have no business being there
change_idx = find(string(sp_lab_c) == 'DWL' &  ismember(min_pia_i, [1,2,4]) &  min_pia_dist <=50);
change_label = vertcat(change_label, horzcat(change_idx, min_pia_i(change_idx)));
%TODO remove duplicates 


labels_v = unique(change_label(:,2));

for x = labels_v' 
    idxs = change_label(find(change_label(:,2) == x),1);
    % x + 1 because the centers don't include none or unmarked and the
    % naming is 0 based
    sp_all(ismember(L, idxs')) = x+1;
    sp_lab_c(idxs) = cellstr(bwlab2layer_map(x+1));
end




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% using the pia dist find misassigned 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% translate from text to bwlabel  
lab_bwlabel = cellfun(@(x) layer2bwlab_map(x), sp_lab_c);

% actual_dist_pia description %%
% 1st col current label
% 2nd col pia_dist from current label
% 3rd col  optimal label by pia dist 
% 4th col pia_dist for optimal
actual_dist_pia = zeros(size(lab_bwlabel,1),3);
for i = 1:size(lab_bwlabel,1)
    lab_col = lab_bwlabel(i)-1;
    if ismember(lab_col, [1:5]) 
        actual_dist_pia(i,1) = pia_layer_dist(lab_col,i);
        [m, mi] = min(pia_layer_dist(:,i));
        actual_dist_pia(i,2) = mi; 
        actual_dist_pia(i,3) = m; 
    else 
        actual_dist_pia(i,:) = NaN;
    end
end
actual_dist_pia = horzcat(lab_bwlabel, actual_dist_pia);


% change the labels as above, do not allow relabeling of the IGL or DWL
% the distance between the pia is too variable 
change_idx = find(actual_dist_pia(:,2) > 40 & ismember(actual_dist_pia(:,1), [1,2,4]));
change_label = horzcat(change_idx, actual_dist_pia(change_idx,3));
% double check the IF distance from center
% get the distance from center to the new label for each considered superpixel 
if_new_d = arrayfun(@(bw_i, label_i) if_dist(label_i, bw_i), change_label(:,1), change_label(:,2));
change_label = change_label(if_new_d <= 0.5,:);

labels_v = unique(change_label(:,2));

for x = labels_v' 
    idxs = change_label(find(change_label(:,2) == x),1);
    % x + 1 because the centers don't include none or unmarked and the
    % naming is 0 based
    sp_all(ismember(L, idxs')) = x+1;
    sp_lab_c(idxs) = cellstr(bwlab2layer_map(x+1));
end   

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Region adjacency filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% take the non-matching elements and check what is next to them 

% start with the existing non matching, sp_all == 1



% identify connected pieces for combination of labels




% oEGL and iEGL

group_all = zeros(size(sp_all,1), size(sp_all,2));
for l = 1:6 
    layer_group = bwlabel(sp_all == l,4);
    
    start_n = max([max(group_all), 1]); 
    group_all = group_all +  (((layer_group > 0)*start_n) + layer_group);     
end 
    
bwlabel(sp_all ==2)


% get regionadjacency 
% get distance to it's nearest neighbor, overall area 



%oEGL_bwlab =  bwlabel(ismember(L, find(string(sp_lab_c) == 'oEGL')),4);
%IGL_bwlab = bwlabel(ismember(L, find(string(sp_lab_c) == 'IGL' )),4);
%ML_bwlab = bwlabel(ismember(L, find(string(sp_lab_c) == 'ML' )),4);
%DWL_bwlab = bwlabel(ismember(L, find(string(sp_lab_c) == 'DWL' )),4);


	

% no gap between the oEGL and the iEGL

% the oEGL shouldn't directly contact the IGL 



% must be at least 1 pixel iEGL between molecular layer






layer =mean_im_cell;



return




end 
