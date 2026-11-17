function layer = cerebellum_superpixels4(file)  

a = imread(file,1);
b = imread(file, 2);
c = imread(file, 3);

a_n = a;
a_n(a_n< 8000) =NaN;

b_n = b;
b_n(b_n< 8000) = NaN;

c_n = c;
c_n(c_n< 6000) = NaN;

b_icor = illumination_cor(b);
% this part also removes grime or bleach marks 
a_nf = adapthisteq(mat2gray(a),'clipLimit',0.01,'Distribution','rayleigh', 'NumTiles', [50 50], 'Range', 'original');
b_nf = adapthisteq(mat2gray(b_icor),'clipLimit',0.01,'Distribution','rayleigh', 'NumTiles', [50 50], 'Range', 'original');
c_nf = adapthisteq(mat2gray(c_n),'clipLimit',0.01,'Distribution','rayleigh', 'NumTiles', [50 50], 'Range', 'original');

a_nfg =  imgaussfilt(a_nf,5);
b_nfg =  imgaussfilt(b_nf,5);
c_nfg =  imgaussfilt(c_nf,5);

col(:,:,1) = a_nf;
col(:,:,2) = b_nf;
col(:,:,3) = c_nf;

% remove low intensity parts of image 
a_nfg_h = mat2gray(a_nfg, [0.4 1]); 
b_nfg_h = mat2gray(b_nfg, [0.4 1]);
c_nfg_h = mat2gray(c_nfg, [0.4 1]);
a_nfg_l = mat2gray(a_nfg, [0 0.4]); 
b_nfg_l = mat2gray(b_nfg, [0 0.4]);
c_nfg_l = mat2gray(c_nfg, [0 0.4]);


bb = mat2gray(b);
all_cerebellum = bb > 0.1;
%all_cerebellum = imclearborder(all_cerebellum, 4);
all_cerebellum = imfill(all_cerebellum, 'holes');


b_norm = (double(b_nfg) - double(a_nfg)) .* all_cerebellum;
b_norm =  mat2gray(b_norm) .* all_cerebellum;
b_norm =  mat2gray(b_norm, [0.4 1]);


dwl_norm = (double(a_nfg_l) - double(c_nfg_l)) ; 
dwl_norm = mat2gray(dwl_norm) .* all_cerebellum;
dwl_norm =  mat2gray(dwl_norm, [0.4 1]);

a_norm = (double(a_nfg_h) - double(c_nfg_h)) ;  
a_norm = mat2gray(a_norm) .* all_cerebellum ;
a_norm =  mat2gray(a_norm, [0.4 1]);


c_norm = (c_nfg_h - (b_norm + a_norm)) ;
c_norm = mat2gray(c_norm) .* all_cerebellum;
c_norm = ordfilt2(c_norm,25,ones(5,5));
c_norm =  imgaussfilt(c_norm,10);
c_norm = ordfilt2(c_norm,25,ones(5,5));
c_norm =  imgaussfilt(c_norm,10);
c_norm = ordfilt2(c_norm,25,ones(5,5));
c_norm =  imgaussfilt(c_norm,10);
c_norm =  mat2gray(c_norm, [0.4 1]);

wml_norm = ((1-b_nfg_l) - c_norm - a_norm - b_norm) ;
wml_norm = mat2gray(wml_norm) .* all_cerebellum ;
wml_norm =  mat2gray(wml_norm, [0.4 1]);

ml_norm = wml_norm - dwl_norm;
ml_norm = mat2gray(ml_norm) .* all_cerebellum ;
ml_norm =  mat2gray(ml_norm, [0.4 1]);
%a_bin = imbinarize(a_nfg, 'adaptive', 'Sensitivity', 0.6);
%c_bin = imbinarize(c_nf, 'adaptive', 'Sensitivity', 0.6);
%a_bin = imbinarize(a_nf, 'adaptive', 'Sensitivity', 0.6);

%red_merge(:,:,1) =  b_norm;
%red_merge(:,:,2) =  c_norm;
%green_merge(:,:,1) =  a_norm;
%green_merge(:,:,2) =  dwl_norm;

%stag_col(:,:,1) = max(red_merge, [], 3);
%stag_col(:,:,2) = max(green_merge, [], 3);
%stag_col(:,:,3) = ml_norm;

cmyk_col(:,:,1) = a_norm;
cmyk_col(:,:,2) = b_norm;
cmyk_col(:,:,3) = wml_norm;
cmyk_col(:,:,4) = c_norm;

%inprof = iccread('USSheetfedCoated.icc');
%outprof = iccread('sRGB.icm');
%C = makecform('icc',inprof,outprof);
C = makecform('cmyk2srgb'); 
I_rgb = applycform(cmyk_col,C);
imshow(I_rgb)


[L,N] = superpixels(I_rgb,20000);
BW = boundarymask(L);
imshow(imoverlay(I_rgb,BW))



%col_g(:,:,1) = a_norm;
%%col_g(:,:,2) = b_norm;
%col_g(:,:,3) = c_norm;


%[L,N] = superpixels(col_g,20000);
sp_vals = regionprops(L,'PixelIdxList','PixelList', 'centroid');




image_cell = {a_norm, b_norm, c_norm, ml_norm, dwl_norm, double(all_cerebellum)}; 


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
c_thresh = V(ceil(size(V,2)*0.25));
%c_sp = find(mean_im(3,:) > c_thresh &  mean_im(2,:) <  0.1);

%ML about 18% 
V = sort(mean_im_trim(4,:), 'descend');
ml_thresh = V(ceil(size(V,2)*0.18));
%ml_sp =  find(mean_im(4,:) > 0.2 & mean_im(1,:) < 0.3 &  mean_im(2,:) < 0.3);

%Deep white layer
V = sort(mean_im_trim(5,:), 'descend');
dwl_thresh = V(ceil(size(V,2)*0.12));


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

% create a mapping from bwlabel to structure 
keySet =   [0,1,2,3,4,5,6];
valueSet = {'none','unmarked', 'iEGL','oEGL','IGL','ML','DWL'};
bwlab2layer_map = containers.Map(keySet,valueSet);
layer2bwlab_map = containers.Map(valueSet,keySet);

% skeletonize the outer egl to get a the outer sufcace 

all_outline = bwmorph(all_cerebellum,'remove');
se = strel('disk',15);
all_outline= imclose(all_outline, se);
all_outline_mask = bwmorph(all_outline,'thicken', 40);
all_outline_mask = bwmorph(all_outline_mask,'bridge', 5);
se = strel('disk',5);
all_outline_mask= imclose(all_outline_mask, se);
all_outline_skel = bwmorph(all_outline,'skel',Inf);

se = strel('disk',2);
egl_in = imbinarize(b_norm+a_norm-c_norm-ml_norm-dwl_norm, 'adaptive', 'Sensitivity', 0.7) .* all_cerebellum;
%egl_in = imfill(egl_in, c_norm > );

% wavelet filter for contour enhancement
Step_Filter_grammes = 25; %size
Step_Filter_Dbw = [4 4]; %vector of width of main wave
Step_Filter_step = pi/12;  %angle step
Step_Filter_DTheta = [0:Step_Filter_step:pi-(Step_Filter_step)]; %angle vector
Step_Filter_type = 0; % type of step filter in {0,1} : 1 step function, 0 Polynomial filter

grammes = Step_Filter_grammes;
Dbw = Step_Filter_Dbw;
DTheta = Step_Filter_DTheta;

[b_wav4]= step_pol_Filtering(b_nf,grammes,DTheta,Dbw,Step_Filter_type);

Step_Filter_Dbw = [3 3]; %vector of width of main wave
Dbw = Step_Filter_Dbw;
[b_wav3]= step_pol_Filtering(b_nf,grammes,DTheta,Dbw,Step_Filter_type);

Step_Filter_Dbw = [5 5]; %vector of width of main wave
Dbw = Step_Filter_Dbw;
[b_wav5]= step_pol_Filtering(b_nf,grammes,DTheta,Dbw,Step_Filter_type);


b_wav(:,:,1) = b_wav3;
b_wav(:,:,2) = b_wav4;
b_wav(:,:,3) = b_wav5;

b_wav_max = max(b_wav,[],3);

b_wav_bin = b_wav_max > 15; 
bwmorph(all_outline,'skel',Inf);



%Step_Filter_Dbw = [6 6]; %vector of width of main wave
%Dbw = Step_Filter_Dbw;
%[b_wav6]= step_pol_Filtering(b_nf,grammes,DTheta,Dbw,Step_Filter_type);



b_wav = mat2gray(b_wav3 + b_wav4 + b_wav4);




% identifying the pial boundary is a little tricky 
% the position has to exact in folia, dividing the outerEGL using the
% skeleton isn't perfect
% using phase symmetry we can enrich the line enough to get a signal

pia = b_norm+a_norm-c_norm-ml_norm-dwl_norm > 0.3;
se = strel('disk',10);
pia = imopen(pia, se);
se = strel('disk',15);
pia = imclose(pia, se);
se = strel('disk',20);
pia = imdilate(pia, se);
pia = imdilate(pia, se);
%se = strel('disk',20);
%pia = imclose(pia, se);

all_comp = b_nf > 0.5 | a_nf > 0.5 ;
%all_comp = all_comp  - (phaseSym > 0);

tt = pia .* all_comp;
tt = pia - tt;
tt = (tt .* all_cerebellum);
se = strel('disk',10);

ml_bin = (ml_norm > 0.6 | c_norm > 0.2);
se = strel('disk',5);
ml_mask = imdilate(ml_bin , se);
se = strel('disk',5);
ml_mask = imopen(ml_mask, se);
se = strel('disk',5);
ml_mask = imdilate(ml_mask , se);

a_bin = (a_norm > 0.4);
se = strel('disk',5);
a_mask = imdilate(a_bin , se);

c_bin = c_norm > 0.2;
%c_mask = imeopen(c_bin, strel('disk',10));
se = strel('disk',20);
c_mask = imdilate(c_bin , se);


tt = tt .* (1-a_mask) .* (1-ml_mask);
%tt = bwmorph(tt,'thicken', 5);
tt = tt | all_outline;
a = bwmorph(tt,'thicken', 5);

imtool(b_norm)


imshow(imcolse(tt,se))
%pia = imbinarize(b_norm+a_norm-c_norm-ml_norm-dwl_norm, 'adaptive', 'Sensitivity', 0.4);


pia = (phaseSym > 0) .* pia;
pia = bwareaopen(pia,10);
pia= imclose(pia, se);



folia_skel = bwmorph(pia,'skel',Inf);
%folia_skel = bwareaopen(folia_skel,200);
folia_skel = folia_skel .* (1-all_outline_mask);

pia_outline = all_outline + folia_skel;
pia_outline = bwmorph(pia_outline,'skel',Inf);

%pia_outline = bwareaopen(pia_outline,500);
    
%pia_outline = filledgegaps(pia_outline, 50);
%pia_outline = bwareaopen(pia_outline,100);
%[edgelist edgeim, etype] = edgelink(pia_outline == 1);
%remove_pia = ismember(edgeim, find(etype == 1));

pia_outline = pia_outline .* (1-remove_pia);
imshow(pia_outline)





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
