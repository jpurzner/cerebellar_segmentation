function layer = cerebellum_superpixels6(file)  

a = imread(file,1);
b = imread(file, 2);
c = imread(file, 3);

%b_icor = illumination_cor(b);
% this part also removes grime or bleach marks 
a_nf = adapthisteq(mat2gray(a),'clipLimit',0.01,'Distribution','rayleigh', 'NumTiles', [50 50], 'Range', 'original');
b_nf = adapthisteq(mat2gray(b),'clipLimit',0.01,'Distribution','rayleigh', 'NumTiles', [50 50], 'Range', 'original');
c_nf = adapthisteq(mat2gray(c),'clipLimit',0.01,'Distribution','rayleigh', 'NumTiles', [50 50], 'Range', 'original');
disp('histogram equilization done')

a_nf = background_subtract(a_nf);
b_nf = background_subtract(b_nf);
c_nf = background_subtract(c_nf);

a_nfg =  imgaussfilt(a_nf,5);
b_nfg =  imgaussfilt(b_nf,5);
c_nfg =  imgaussfilt(c_nf,5);
disp('gaussian filtering done')

% anisotropic diffusion for phasesym 
b_nfad = anisodiff(b_nf,20, 20, 0.25, 1);
a_nfad = anisodiff(a_nf,20, 20, 0.25, 1);
disp('anisotropic diffusion done')

% phasesym for iden
[b_phasesym, orientation, totalEnergy, T] = phasesym(b_nfad, 5, 6, 3, 2.5, 0.55, 2.0, 0);
%[b_ps_nmxs,location] = nonmaxsup(b_phasesym, orientation, 3);

[a_phasesym, orientation, totalEnergy, T] = phasesym(a_nfad, 5, 6, 3, 2.5, 0.55, 2.0, 0);
%[a_ps_nmxs,location] = nonmaxsup(a_phasesym, orientation, 3);
disp('phase symmetry done')

col(:,:,1) = a_nf;
col(:,:,2) = b_nf;
col(:,:,3) = c_nf;

% remove low intensity parts of image 
a_nfg_h = background_subtract(mat2gray(a_nfg, [0.4 1])); 
b_nfg_h = background_subtract(mat2gray(b_nfg, [0.4 1]));
c_nfg_h = background_subtract(mat2gray(c_nfg, [0.4 1]));
a_nfg_l = mat2gray(a_nfg, [0 0.4]); 
b_nfg_l = mat2gray(b_nfg, [0 0.4]);
c_nfg_l = mat2gray(c_nfg, [0 0.4]);



all_cerebellum = b_nfad  + a_nfad;
all_cerebellum = background_subtract(all_cerebellum);
all_cerebellum = all_cerebellum > 0.4; 
cerebellum_boundary = bwboundaries(all_cerebellum, 'noholes');
[~,b_i] = sort(cellfun(@length,cerebellum_boundary), 'descend');
all_outline = zeros(size(b_nfad));
boundary_y = cerebellum_boundary{b_i(1)}(:,1);
boundary_x = cerebellum_boundary{b_i(1)}(:,2);
% smooth here 
for n = 1:size(boundary_y,1) 
    all_outline(boundary_y(n), boundary_x(n)) = 1;
end

all_cerebellum = imclose(all_outline, strel('disk', 10));
all_cerebellum = imfill(all_cerebellum, 'holes');
all_outline = bwperim(all_cerebellum);
all_outline_mask = imdilate(all_outline,strel('disk', 50));

b_norm = (double(b_nfg) - double(a_nfg)) .* all_cerebellum;
b_norm =  mat2gray(b_norm) .* all_cerebellum;
b_norm =  background_subtract(b_norm);

a_norm = (double(a_nfg_h) - (double(c_nfg_h) .* (1 - all_outline_mask))) ;  
a_norm = mat2gray(a_norm) .* all_cerebellum ;
a_norm =  background_subtract(a_norm);

c_norm = (c_nfg_h - (a_norm + b_norm)) ;
c_norm = mat2gray(c_norm) .* all_cerebellum;
c_norm = ordfilt2(c_norm,25,ones(5,5));
c_norm =  imgaussfilt(c_norm,10);
c_norm = ordfilt2(c_norm,25,ones(5,5));
c_norm =  imgaussfilt(c_norm,10);
c_norm = ordfilt2(c_norm,25,ones(5,5));
c_norm =  imgaussfilt(c_norm,10);
c_norm = background_subtract(c_norm);

% take 30% of cerebellum 
c_bin = thresh_by_area( c_norm, all_cerebellum, 0.3 );
c_bin = c_bin .* (1 - all_outline_mask);
c_bin = bwareafilt(c_bin > 0,[5000 Inf]);
c_mask = imclose(c_bin, strel('disk',50));
c_mask = imdilate(c_mask, strel('disk',8));


% revise a and b norm with c norm 
a_norm = a_norm - c_norm;  
a_norm(a_norm < 0) = 0; 
a_bin = thresh_by_area( a_norm, all_cerebellum, 0.1 );

b_norm = b_norm - c_norm;  
b_norm(b_norm < 0) = 0;
b_bin = thresh_by_area( b_norm, all_cerebellum, 0.1 );


dwl_norm = (double(a_nfg_l) - double(c_nfg_l) - double(a_nfg_h) - (double(a_nfg_h)))  ; 
dwl_norm = mat2gray(dwl_norm) .* all_cerebellum;
dwl_norm = ordfilt2(dwl_norm,25,ones(5,5));
dwl_norm =  imgaussfilt(dwl_norm,5);
dwl_norm = background_subtract(dwl_norm);


dwl_bin = thresh_by_area( dwl_norm, all_cerebellum, 0.2 );
dwl_bin = dwl_bin .* (1 - imdilate(all_outline_mask,strel('disk', 70)));
dwl_bin = dwl_bin .* (1 - c_bin);
dwl_mask = imdilate(dwl_bin, strel('disk',2));
dwl_mask = imclose(dwl_mask, strel('disk',15));
dwl_mask = imdilate(dwl_mask, strel('disk',10));
dwl_mask = bwareafilt(dwl_mask > 0,[5000 Inf]);
dwl_bwl = bwlabel(dwl_mask);
dwl_bwl_counts = histcounts(dwl_bwl, max(unique(dwl_bwl)));
[~, dwl_max_lab] = max(dwl_bwl_counts(2:end));
dwl_mask = dwl_bwl == dwl_max_lab;


wml_norm = ((1-b_nfg_l) - c_norm - a_norm) ;
wml_norm = mat2gray(wml_norm) .* all_cerebellum ;
wml_norm = background_subtract(wml_norm);

ml_norm = (wml_norm - dwl_norm) .* (1 - imdilate(all_outline,strel('disk', 15)));
ml_norm = mat2gray(ml_norm) .* all_cerebellum ;
ml_norm = background_subtract(ml_norm);
ml_bin = thresh_by_area( ml_norm, all_cerebellum, 0.2 );

montage([a_norm b_norm c_norm; ...
    ml_norm dwl_norm all_outline;...
    a_bin b_bin c_bin; ...
    ml_bin dwl_bin all_cerebellum])

%% -----------------------------------------------------------------------
% FIND THE PIAL BOUNDARY 
% identifying the pial boundary is a little tricky 
% the position has to exact in folia, dividing the outerEGL using the
% skeleton isn't perfect
% using phase symmetry we can enrich the line enough to get a signal



% use EGL channels to define the approximate location of the pia within the
% buried portion of the folia 
%pia = imbinarize(mat2gray(b_norm+a_norm-c_norm-ml_norm-dwl_norm), 'adaptive', 'Sensitivity',0.8);

pia_norm = b_norm+a_norm-c_norm-ml_norm -dwl_norm;
pia_norm(pia_norm < 0) = 0; 
pia_bin = thresh_by_area( pia_norm, all_cerebellum, 0.12 );

pia = imdilate(pia_bin, strel('disk',10));
pia = imerode(pia, strel('disk',5));;
pia = bwareafilt(pia > 0,[20000 Inf]);
pia = imdilate(pia, strel('disk',10));
pia_thin = bwmorph(pia,'thin', Inf);


se = strel('disk',15);
pia_mask = imdilate(pia_thin, se);

pia_phasesym = ((b_phasesym + a_phasesym) .* pia_mask);
[pia_line,location] = nonmaxsup(pia_phasesym, orientation, 3);
pia_line = pia_line > 0.5; 
pia_line = bwareafilt(pia_line > 0,[40 Inf]);
%pia_bin2 = pia_phasesym >0.5 ;

pia_inner = pia_line .* (1-all_outline_mask);
%
csf_space = ((1- b_nfad)  > 0.5) .* all_cerebellum;
pia_csf = (csf_space + pia_inner) .* (1-all_outline_mask); 
pia_csf_lab = bwlabel(pia_csf);

%filter by area of overlap 
all_pia_ov = unique(pia_csf_lab .* pia_inner);
all_pia_ov = sort(all_pia_ov);
all_pia_ov = all_pia_ov(2:end);

pia_ov_t = histc(reshape(pia_csf_lab .* pia_inner,[],1), all_pia_ov) ./ histc(reshape(pia_csf_lab ,[],1), all_pia_ov);

pia_ov_sel = all_pia_ov(pia_ov_t > 0.001);

pia_deep = (1 - (a_nfg > 0.4 )) .* ismember(pia_csf_lab, pia_ov_sel);
pia_deep = pia_deep | pia_inner;
pia_deep = bwareafilt(pia_deep > 0,[20 Inf]);

pia_deep_thin = bwmorph(pia_deep .* (1- c_mask),'thin', Inf);
pia_deep_thin = filledgegaps(pia_deep_thin, 101);

pia_deep = imdilate(pia_deep | pia_deep_thin,strel('disk',3));
pia_deep = imerode(pia_deep,strel('disk',2));
pia_deep = imfill(pia_deep, 'holes');
tt = bwmorph(pia_deep,'thin', Inf);
tt = filledgegaps(tt, 10);

pia_deep = bwareafilt(pia_deep > 0,[500 Inf]);

pia_outline = all_outline | pia_deep;
%pia_outline_thin = bwmorph(pia_outline,'thin', Inf);

% for every pia_deep_segment get the start and end and find the min
% distance to the the pia outline then burn a line to this point 
pia_deep_bwl = bwlabel(pia_deep);  

% the dwl couldn't be subtracted initially because it distorted images 

pia_discard = unique(pia_deep_bwl .* c_mask);
[outline_y, outline_x] = find(all_outline);  % x and y are column vectors.
outline_xy = [outline_x, outline_y];

pia_junction = zeros(size(pia_outline)); 
for n = 1:max(unique(pia_deep_bwl))
    if ~ismember(n,pia_discard)
        [~, ~, re, ce] = findendsjunctions(pia_deep_bwl == n);
        deep_ends = [re, ce];
        [D,I] = pdist2(outline_xy, deep_ends, 'euclidean', 'Smallest',1);
        [min_D, deep_end_I] = min(D);
        if min_D < 200 
            outline_closest = outline_xy(I(deep_end_I),:);
            deep_ends_closest = deep_ends(deep_end_I,:);
            conn_line = [outline_closest(2), outline_closest(1), deep_ends_closest(2), deep_ends_closest(1)];
            pia_junction = insertShape(double(pia_junction),'line',conn_line);
        end
    end
end 


pia_deep(ismember(pia_deep_bwl, pia_discard)) = 0;
pia_outline = all_outline | pia_deep;
pia_outline = pia_outline | (pia_junction(:,:,1) > 0);
pia_outline = bwareafilt(pia_outline > 0,[100 Inf]);

cmyk_col(:,:,1) = a_norm .* (1 - pia_outline) .* all_cerebellum;
cmyk_col(:,:,2) = b_norm .* (1 - pia_outline) .* all_cerebellum;
cmyk_col(:,:,3) = wml_norm .* (1 - pia_outline) .* all_cerebellum;
cmyk_col(:,:,4) = c_norm .* (1 - pia_outline) .* all_cerebellum;

%inprof = iccread('USSheetfedCoated.icc');
%outprof = iccread('sRGB.icm');
%C = makecform('icc',inprof,outprof);
C = makecform('cmyk2srgb'); 
I_rgb = applycform(cmyk_col,C);

figure
imshow(I_rgb)


[L,N] = superpixels(I_rgb,20000);
BW = boundarymask(L);
imshow(imoverlay(I_rgb,BW))



%col_g(:,:,1) = a_norm;
%%col_g(:,:,2) = b_norm;
%col_g(:,:,3) = c_norm;


%[L,N] = superpixels(col_g,20000);
sp_vals = regionprops(L,'PixelIdxList','PixelList', 'centroid');




image_cell = {a_norm, b_norm, c_norm, ml_norm, dwl_norm, double(all_cerebellum), double(pia_outline)}; 


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


%M = containers.Map(1:length(k),k, 'UniformValues',false, 'KeyType', 'uint16', 'ValueType', 'uint16');
%M = containers.Map('KeyType', 'uint32', 'ValueType', 'uint16');
%M = containers.Map(uint32(1:length(k)),uint16(k), 'UniformValues',true);

% start with heuristics to get the best matching superpixels then fill in
% remainder by proximity

% the thresholds still vary from sample to sample
% start by taking the top x% of superpixels and assigning 
mean_im_trim = mean_im(:,mean_im(6,:) ==1);

%iGL about 8% of trim 

scatter(log2(mean_im_trim(1,:)./mean_im_trim(2,:)),  log2(mean_im_trim(3,:)./mean_im_trim(4,:)) )
scatter(mean_im_trim(1,:),mean_im_trim(2,:))

%outer EGL
V = sort(mean_im_trim(2,:), 'descend');
b_thresh = V(ceil(size(V,2)*0.10));

%inner EGL about 8% of trim 
V = sort(mean_im_trim(1,:), 'descend');
a_thresh = V(ceil(size(V,2)*0.10));
%a_sp = find(mean_im(1,:) >= thresh &  mean_im(2,:) < 0.2);

%IGL about 25% of trim 
V = sort(mean_im_trim(3,:), 'descend');
c_thresh = V(ceil(size(V,2)*0.30));
%c_sp = find(mean_im(3,:) > c_thresh &  mean_im(2,:) <  0.1);

%ML about 18% 
V = sort(mean_im_trim(4,:), 'descend');
ml_thresh = V(ceil(size(V,2)*0.18));
%ml_sp =  find(mean_im(4,:) > 0.2 & mean_im(1,:) < 0.3 &  mean_im(2,:) < 0.3);

%Deep white layer , artificially supress this to prevent takeover of the ML
% IGL junction 
V = sort(mean_im_trim(5,:), 'descend');
dwl_thresh = V(ceil(size(V,2)*0.06));


a_sp = find(mean_im(1,:) >= a_thresh & mean_im(2,:) < b_thresh ... 
                                     & mean_im(3,:) < c_thresh ...
                                     & mean_im(4,:) < ml_thresh ...
                                     & mean_im(5,:) < dwl_thresh);

b_sp = find(mean_im(2,:) >= b_thresh & mean_im(1,:) < a_thresh ... 
                                     & mean_im(3,:) < c_thresh ...
                                     & mean_im(4,:) < ml_thresh ...
                                     & mean_im(5,:) < dwl_thresh);

c_sp = find(mean_im(3,:) >= c_thresh & mean_im(1,:) < a_thresh ... 
                                     & mean_im(2,:) < b_thresh ...
                                     & mean_im(4,:) < ml_thresh ...
                                     & mean_im(5,:) < dwl_thresh);
                                 
ml_sp =  find(mean_im(4,:) >= ml_thresh  & mean_im(1,:) < a_thresh ... 
                                         & mean_im(2,:) < b_thresh ...
                                         & mean_im(3,:) < c_thresh ...
                                         & mean_im(5,:) < dwl_thresh);
                                 
dwl_sp =  find(mean_im(5,:) > dwl_thresh & mean_im(1,:) < a_thresh ... 
                                         & mean_im(2,:) < b_thresh ...
                                         & mean_im(3,:) < c_thresh ...
                                         & mean_im(4,:) < ml_thresh);

                                     
                                     
                                     
                                     
pia_sp = setdiff(find(mean_im(7,:) == 1), horzcat(a_sp,b_sp));                                     
                                     
% cerebellum mask in superpixels 
cbl_sp = find(mean_im(6,:) == 1);



%re-lable the superpixels 
sp_all = zeros(size(L,1),size(L,2));
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
sp_all(ismember(L, pia_sp)) = 7;
sp_lab_c(pia_sp) = cellstr('pia');

% create a mapping from bwlabel to structure 
keySet =   [0,1,2,3,4,5,6,7];
valueSet = {'none','unmarked', 'iEGL','oEGL','IGL','ML','DWL','pia'};
bwlab2layer_map = containers.Map(keySet,valueSet);
layer2bwlab_map = containers.Map(valueSet,keySet);

sp_all(pia_outline == 1) = 7;

figure
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
change_idx = find(string(sp_lab_c) == 'unmarked' &  min_pia_dist <=30);
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
% RETURN results

imagesc(sp_all)

layer =mean_im_cell;
return





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

%%%%%%%%% not used 
% change the labels as above, do not allow relabeling of the IGL or DWL
% the distance between the pia is too variable 
%change_idx = find(actual_dist_pia(:,2) > 40 & ismember(actual_dist_pia(:,1), [1,2,4]));
%%%%%%%%


%allow reassignment of IGL and IGL 
change_idx = find(actual_dist_pia(:,2) > 40);  %& ismember(actual_dist_pia(:,1), [1,2,4]));
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











end 
