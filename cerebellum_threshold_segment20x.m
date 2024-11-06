function layer = cerebellum_threshold_segment20x(file, varargin)  

%segment the developing cerebellum at P7 into multiple layers from 
% the outside in.  The layer are as follows  :
%    oEGL: DAPI only outer layer  
%    iEGL: p27 
%    ML: 
% 

% file must contain images with DAPI, p27 and NeuN 
% pixel_size: size of pixel in micometers, at 20x was 0.5119049, at 10x was
%               1.0239
% scale_factor: the constant that determines the size of the image```
%               operations, 1 for 10x 2 for 20x, if not provided then 
%               determined by pixel_size if provided and will default 
%               to 1 if pixel_size is not given
% TODO specify the order of the images 



% Set up default parameter values as needed
if nargin < 2
    pixel_size  = 0.5119049;
    %pixel_width <- 0.5119049  # Width of a pixel in micrometers
end
% optionally overide the scalling factor 
if nargin < 3
	scaling = 1;
end
    

% the image order for these files 

a = imread(file,1);
c = imread(file, 2);
b = imread(file, 3);

%a = illumination_cor(a);
%b = illumination_cor(b);
%c = illumination_cor(c);
% this part also removes grime or bleach marks 
a_nf = adapthisteq(mat2gray(a),'clipLimit',0.01,'Distribution','rayleigh', 'NumTiles', [50 50], 'Range', 'original');
b_nf = adapthisteq(mat2gray(b),'clipLimit',0.01,'Distribution','rayleigh', 'NumTiles', [50 50], 'Range', 'original');
c_nf = adapthisteq(mat2gray(c),'clipLimit',0.01,'Distribution','rayleigh', 'NumTiles', [50 50], 'Range', 'original');
disp('histogram equilization done')

% keep for later
a_set_thresh = a_nf ; 
b_set_thresh = b_nf ; 

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
[a_ps_nmxs,location] = nonmaxsup(a_phasesym, orientation, 3);
disp('phase symmetry done')

a_level = graythresh(a_nf);
b_level = graythresh(b_nf);
c_level = graythresh(c_nf);

col(:,:,1) = a_nf;
col(:,:,3) = b_nf;
col(:,:,2) = c_nf;

% remove low intensity parts of image 
a_nfg_h = background_subtract(mat2gray(a_nfg, [a_level 1])); 
b_nfg_h = background_subtract(mat2gray(b_nfg, [b_level 1]));
c_nfg_h = background_subtract(mat2gray(c_nfg, [c_level 1]));
a_nfg_l = mat2gray(a_nfg, [0 a_level]); 
b_nfg_l = mat2gray(b_nfg, [0 b_level]);
c_nfg_l = mat2gray(c_nfg, [0 c_level]);

%-------------------------------------------------------------------------
% the Purkinje cells are creating a iEGL like layer in the final image
a_pc = a_nf;
a_pc(a_nf < 0.25 | a_nf > 0.8) = 0;
a_pc = mat2gray(a_pc);

pk = (a_nf - b_nf).*(a_nf -c_nf) .* a_nf; 
pk = adapthisteq(pk);
%pk = pk + mat2gray(a_pc);
pc_bin  = pk > 0.1;
pc_bin = bwareafilt(pc_bin > 0,[40 Inf]);
pc_bin = imfill(pc_bin,'holes');
pc_bwl = bwlabel(pc_bin); 
% code snippet from image analyst 
% https://www.mathworks.com/matlabcentral/answers/111855-how-to-use-the-imfilter-function-to-pick-up-blobs-that-are-circular-linear
measurements = regionprops(pc_bwl, 'Area', 'Perimeter');
allAreas = [measurements.Area];
allPerimeters = [measurements.Perimeter];
circularities = allPerimeters .^ 2 ./ (4 * pi * allAreas);
roundObjectsIndexes = find(circularities < 3 & allAreas > 40 & allAreas < 300 ); % or whatever value you want.
% Extract only those blobs that meet our criteria, and
% eliminate those blobs that don't meet our criteria.
% Note how we use ismember() to do this.
keeperBlobsImage = ismember(pc_bwl, roundObjectsIndexes) ; % A labeled image
pc_bin = keeperBlobsImage > 0;
pc_bin = imopen(pc_bin, strel('disk',2 * scaling));
pc_bin = bwareafilt(pc_bin > 0,[40 Inf]);
pc_bin = imdilate(pc_bin,strel('disk', 4 * scaling));

%imshow(label2rgb(keeperBlobsImage, 'winter', 'k'))
%imshow(imoverlay(a, bwperim(  pc_bin)))

clear measurements allAreas allPerimeters circularities roundObjectsIndexes
clear keeperBlobsImage
%clear a_nf b_nf c_nf a b c 

all_cerebellum = b_nfad  + a_nfad;
all_cerebellum = background_subtract(all_cerebellum);
all_cerebellum = all_cerebellum > 0.1; 
% handles the case when there is a breach in the EGL 
all_cerebellum_large = imclose(all_cerebellum, strel('disk', 20*scaling));
all_cerebellum_large = imfill(all_cerebellum_large, 'holes');
all_cerebellum_large = all_cerebellum_large - all_cerebellum;
all_cerebellum_large = imclose(all_cerebellum_large, strel('disk', 5*scaling));
all_cerebellum_large = imfill(all_cerebellum_large, 'holes');
all_cerebellum_large = bwareafilt(all_cerebellum_large > 0,[5000 Inf]);
all_cerebellum = all_cerebellum | all_cerebellum_large;

cerebellum_boundary = bwboundaries(all_cerebellum, 'noholes');
[~,b_i] = sort(cellfun(@length,cerebellum_boundary), 'descend');
all_outline = zeros(size(b_nfad));
boundary_y = cerebellum_boundary{b_i(1)}(:,1);
boundary_x = cerebellum_boundary{b_i(1)}(:,2);
% smooth here 
for n = 1:size(boundary_y,1) 
    all_outline(boundary_y(n), boundary_x(n)) = 1;
end
all_cerebellum = imclose(all_outline, strel('disk', 10*scaling));
all_cerebellum = imfill(all_cerebellum, 'holes');
all_outline = bwperim(all_cerebellum);
all_outline_mask = imdilate(all_outline,strel('disk', 50*scaling));
disp('cerebellum outline done')


c_norm = c_nfg_h .* (1- all_outline_mask) ;
c_norm = mat2gray(c_norm) .* all_cerebellum;
c_norm = ordfilt2(c_norm,25,ones(5,5));
c_norm =  imgaussfilt(c_norm,10);
c_norm = ordfilt2(c_norm,25,ones(5,5));
c_norm =  imgaussfilt(c_norm,10);
c_norm = ordfilt2(c_norm,25,ones(5,5));
c_norm =  imgaussfilt(c_norm,10);
c_norm = background_subtract(c_norm);

% take 30% of cerebellum 
c_bin = thresh_by_area( c_norm, all_cerebellum, 0.4 );
c_bin = c_bin .* (1 - all_outline_mask);
c_bin = bwareafilt(c_bin > 0,[5000 Inf]);
c_bin = imclose(c_bin, strel('disk',10*scaling));
c_bin = bwareafilt(c_bin > 0,[20000 Inf]);



c_whole = im_get_overlap(c_bin, c_nfg_h > c_level, 0.6, 200);
ca_whole = im_get_overlap(c_bin, a_nfg_h > a_level, 0.6, 200);
cb_whole = im_get_overlap(c_bin, b_nfg_h > b_level, 0.3, 200);
c_final = (ca_whole + c_whole + cb_whole)  > 0 ; 
c_final = imclose(c_final, strel('disk', 5*scaling));
c_final  = smallfill( c_final, 10000 );
c_final = bwareafilt(c_final > 0,[200000 Inf]);


disp('igl done')

%% get the white matter areas deep to the igl 
% they have relatively high levels of p27 comapred to the ML 
dwl_norm = (double(a_nfg_l) - double(c_nfg_l) - double(a_nfg_h) - (double(a_nfg_h)))  ; 
dwl_norm = mat2gray(dwl_norm) .* all_cerebellum;
dwl_norm = ordfilt2(dwl_norm,25,ones(5,5));
dwl_norm =  imgaussfilt(dwl_norm,5);
dwl_norm = background_subtract(dwl_norm);

dwl_bin = thresh_by_area( dwl_norm, all_cerebellum, 0.08 );
dwl_bin = dwl_bin .* (1 - imdilate(all_outline_mask,strel('disk', 70*scaling)));
dwl_bin = dwl_bin .* (1 - c_bin);

dwl_mask = imdilate(dwl_bin, strel('disk',2*scaling));
dwl_mask = imclose(dwl_mask, strel('disk',15*scaling));
dwl_mask = imdilate(dwl_mask, strel('disk',10*scaling));
dwl_mask = bwareafilt(dwl_mask > 0,[4000 Inf]);

igl_in =  c_final | dwl_bin ;
dwl_phase = (a_phasesym .* igl_in > 0.4);
igl_in = igl_in | dwl_phase ;

igl_in = imdilate(igl_in, strel('disk',5*scaling));
igl_in = imclose(igl_in, strel('disk',5*scaling));
igl_in = smallfill( igl_in, 1000 );
igl_in = bwareafilt(igl_in > 0,[70000 Inf]);

% remove stuff outside of egl_mask 
ml_out = (1-igl_in) .* all_cerebellum ;
ml_out = ml_out | all_outline_mask;
ml_out = imdilate(ml_out, strel('disk', 10*scaling)); 
ml_out = imerode(ml_out, strel('disk', 5*scaling)); 
ml_out = bwareafilt(ml_out > 0,[30000 Inf]);
ml_out_deep = ml_out .* (1 - imdilate(all_outline_mask, strel('disk', 60*scaling))) ; 

ml_out_deep = imfill(ml_out_deep, 'holes');
ml_out = ml_out | ml_out_deep;
c_out = im_get_overlap(ml_out, igl_in, 0.5, 100 );
c_out = c_out > 6; 

igl_in = igl_in .* (1 - c_out);
igl_in = igl_in == 1; 
c_final = c_final .* (1 -  c_out);
c_final = c_final == 1;

% TO DO, needs to be more robust to pieces broken off DWL
dwl_final = igl_in - c_final;
dwl_final = dwl_final | dwl_phase;
dwl_final = smallfill( dwl_final, 5000 );
dwl_final = imerode(dwl_final, strel('disk',5));
dwl_final = imdilate(dwl_final, strel('disk',10));
dwl_final = dwl_final .* (1-c_final);
dwl_final = bwareafilt(dwl_final > 0,[20000 Inf]);


igl_in = c_final | dwl_final;
igl_in_c = imclose(igl_in, strel('disk',20));
igl_in = igl_in | igl_in_c;
%igl_in = imfill(igl_in, 'holes');
igl_in  = smallfill( igl_in, 200000 );
dwl_final = igl_in - c_final;


% remove stuff outside of egl_mask 
ml_out = (1-igl_in) .* all_cerebellum ;
ml_out = ml_out | all_outline_mask;
ml_out = imdilate(ml_out, strel('disk', 10*scaling)); 
ml_out = imerode(ml_out, strel('disk', 5*scaling)); 
ml_out = bwareafilt(ml_out > 0,[30000 Inf]);
ml_out_deep = ml_out .* (1 - imdilate(all_outline_mask, strel('disk', 60))) ; 

ml_out_deep = imfill(ml_out_deep, 'holes');
ml_out = ml_out | ml_out_deep;
c_out = im_get_overlap(ml_out, igl_in, 0.5, 100 );
c_out = c_out > 6;



%imshow(imoverlay(dwl_final, c_final, [0 1 0]))

igl_in = igl_in .* (1 - c_out);
igl_in = igl_in == 1; 
c_final = c_final .* (1 -  c_out);
c_final = c_final == 1;
dwl_final = igl_in - c_final;
ml_out = imdilate(ml_out, strel('disk', 3*scaling)); 
dwl_final = dwl_final .* (1 - ml_out); 
dwl_final = bwareafilt(dwl_final > 0,[20000 Inf]);

igl_in = c_final | dwl_final;
dwl_final = igl_in - c_final;

igl_in_mask = imdilate(igl_in, strel('disk',30*scaling));


disp('DWL done')


%% Robustly call the EGL 
%egl_norm = (double(b_nfad) + double(a_nfad)) .* (1- igl_in) .* all_cerebellum;
egl_norm = double(b_nfg)  .* (1- igl_in_mask) .* all_cerebellum;
egl_norm = ordfilt2(egl_norm,25,ones(5,5));
%egl_norm = background_subtract(mat2gray(egl_norm));
egl_bin = thresh_by_area( egl_norm, all_cerebellum, 0.15 );
egl_bin = bwareafilt(egl_bin > 0,[50 Inf]);

egl_strict = thresh_by_area( egl_norm, all_cerebellum, 0.1 );
egl_strict = imdilate(egl_strict, strel('disk',10*scaling));
egl_strict = bwareafilt(egl_strict > 0,[5000 Inf]);



egl_mask = imdilate(egl_bin, strel('disk',10*scaling));
egl_mask = bwareafilt(egl_mask > 0,[5000 Inf]);

egl_thin = bwmorph(egl_mask,'thin', Inf);
[~, egl_thin, ~] = edgelink(egl_thin);
egl_thin = filledgegaps(egl_thin, 101);
egl_thin = bwareafilt(egl_thin > 0,[100 Inf]);
egl_thin = filledgegaps(egl_thin, 201);
egl_thin = imdilate(egl_thin, strel('disk',3*scaling));

egl_mask = egl_mask | egl_thin;
egl_mask = bwareafilt(egl_mask > 0,[20000 Inf]);
egl_mask = imdilate(egl_mask, strel('disk',5*scaling));

egl_phase = im_get_overlap(egl_mask, a_phasesym > 0.1, 0.1, 200 );
egl_mask = imdilate(egl_mask | egl_phase, strel('disk',5*scaling));

%egl_thin = bwmorph(egl_mask,'thin', Inf);
%egl_thin_mask = imdilate(egl_thin, strel('disk',10));

%% refine the purkinje cell zone, which includes more then PC 
pc_layer_mask = igl_in_mask - igl_in - egl_mask;
pc_layer_mask = pc_layer_mask == 1;
pc_layer_bin = (pc_layer_mask .* b_nfg) > 0.4;


a_set_thresh_mean = mean(mean(a_set_thresh));
a_set_thresh_egl = a_set_thresh;
a_set_thresh_egl((1-all_cerebellum) > 0) = a_set_thresh_mean;
%aS_25 = get_thresh_by_area_adapt(a_set_thresh_egl, c_final, 0.2);
%aS_50 = get_thresh_by_area_adapt(a_set_thresh_egl, c_final, 0.35);


T = adaptthresh(mat2gray(a_set_thresh_egl), 0.01,'NeighborhoodSize' , 2*floor(size(a_nfg)/25)+1);
a_bin_1 = imbinarize(double(a_set_thresh),T);




T = adaptthresh(mat2gray(a_set_thresh_egl), 0.1,'NeighborhoodSize' , 2*floor(size(a_nfg)/25)+1);
a_bin_25 = imbinarize(double(a_set_thresh),T);

%a_nf_flat = adapthisteq(mat2gray(a),'clipLimit',0.05,'Distribution','rayleigh', 'NumTiles', [50 50], 'Range', 'original');
%a_set_thresh_mean = mean(mean(a_nf_flat));
%a_set_thresh_egl = a_set_thresh;
%a_set_thresh_egl((1-all_cerebellum) > 0) = a_set_thresh_mean;


T = adaptthresh(mat2gray(a_set_thresh_egl), 0.3,'NeighborhoodSize' , 2*floor(size(a_nfg)/25)+1);
a_bin = imbinarize(double(a_set_thresh),T);
a_bin = bwareafilt(a_bin > 0,[10 Inf]);
a_bin_g = imgaussfilt(double(a_bin),5);




a_bin_smooth = a_bin_g >= 0.2;
a_bin_smooth = a_bin_smooth .* egl_mask;
a_bin_smooth_c = imclose(a_bin_smooth, strel('disk',3*scaling));
a_bin_smooth = smallfill(a_bin_smooth | a_bin_smooth_c, 200);
a_bin_smooth = bwareafilt(a_bin_smooth > 0,[200 Inf]);

%imshow(imoverlay(a_nf, bwperim(a_bin_smooth ),  [0 1 0]))

a_mask = imdilate(a_bin_smooth, strel('disk', 5*scaling));
a_thin = bwmorph(a_mask,'thin', Inf);
a_bin_25 = a_bin_25 | imdilate(a_thin, strel('disk', 1*scaling)); 
a_bin_25 = imgaussfilt(double(a_bin_25),5);
a_bin= a_bin_25 > 0.2;
a_bin= bwareafilt(a_bin_smooth > 0,[500 Inf]);

[~, a_thin, ~] = edgelink(a_thin);
a_thin = filledgegaps(a_thin, 51);
a_thin = filledgegaps(a_thin, 101);
a_thin = bwareafilt(a_thin > 0,[100 Inf]);
a_thin = filledgegaps(a_thin, 151);

imshow(imoverlay(a_nf, bwperim(a_bin ),  [0 1 0]))

%oegl_out = imdilate(a_thin, strel('disk', 5)) ...
%    | imdilate(all_outline, strel('disk',20)) ...
%    | a_bin ...
%    | imdilate(egl_thin, strel('disk',5));


oegl_out = imdilate(a_thin, strel('disk', 5*scaling)) ...
    | imdilate(all_outline, strel('disk',20*scaling)) ...
    | a_bin;
oegl_out = imdilate(oegl_out, strel('disk',5*scaling));
oegl_out = smallfill( oegl_out, 1000000 );
oegl_out = imerode(oegl_out, strel('disk',5*scaling));


b_set_thresh_mean = mean(mean(b_set_thresh));
b_set_thresh_egl = b_set_thresh;
b_set_thresh_egl((1-all_cerebellum) > 0) = b_set_thresh_mean;
%bS_25 = get_thresh_by_area_adapt(b_set_thresh_egl, c_final, 0.2);
%bS_50 = get_thresh_by_area_adapt(b_set_thresh_egl, c_final, 0.35);
T = adaptthresh(mat2gray(b_set_thresh_egl), 0.2,'NeighborhoodSize' , 2*floor(size(a_nfg)/25)+1);
b_bin = imbinarize(double(b_set_thresh),T);
T = adaptthresh(mat2gray(b_set_thresh_egl), 0.01,'NeighborhoodSize' , 2*floor(size(a_nfg)/25)+1);
b_bin_25 = imbinarize(double(b_set_thresh),T);
b_bin_g = imgaussfilt(double(b_bin),5);
b_bin_smooth = b_bin_g >= 0.2;

b_bin_smooth = b_bin_smooth .* egl_mask;
b_bin_smooth = imclose(b_bin_smooth, strel('disk',2));
b_bin_smooth = smallfill(b_bin_smooth | b_bin, 200);
b_bin_smooth = bwareafilt(b_bin_smooth > 0,[100 Inf]);
b_bin_whole = b_bin;
b_bin = b_bin_smooth;
b_bin = b_bin .* oegl_out;


b_norm = b_nf - a_nf - igl_in;
b_norm( b_norm <= 0) = 0;
b_norm = mat2gray(b_norm );
b_bin_steal = thresh_by_area( b_norm, all_cerebellum, 0.05 );
%b_bin = imclose(b_bin,  strel('disk', 2));
%b_bin = bwareafilt(b_bin > 0,[100 Inf]);
b_bin_steal = b_bin_steal .* oegl_out;




%imshow(imoverlay(b_nf, bwperim(b_bin ),  [0 1 0]))
%imshow(imoverlay(b_nf, bwperim(egl_mask ),  [0 1 0]))
%T = adaptthresh(mat2gray(a_set_thresh_egl), 0.1,'NeighborhoodSize' , 2*floor(size(a_nfg)/25)+1);
%a_bin_01 = imbinarize(double(a_set_thresh),T);
%T = adaptthresh(mat2gray(b_set_thresh_egl), 0.1,'NeighborhoodSize' , 2*floor(size(a_nfg)/25)+1);
%b_bin_01 = imbinarize(double(b_set_thresh),T);
%dual_egl_mask = a_bin_01 | b_bin_01;
%egl_vals = dual_egl_mask .* oegl_out .* log2( b_set_thresh ./ a_set_thresh  );
%egl_vals = dual_egl_mask .* oegl_out .* log2( b_nf ./ a_nf  );
%imshow(egl_vals > 0.1)
%imshow(egl_vals > 0.1 )



b_bin = b_bin .* (1- a_bin); 
b_bin = b_bin | b_bin_steal;
b_mask = imclose(b_bin, strel('disk', 10*scaling));
b_mask = imdilate(b_bin, strel('disk', 20*scaling));
b_mask = (b_mask .* (1- a_bin)) | b_bin_steal ; 


ml_layer_mask = imdilate(igl_in_mask, strel('disk', 10*scaling)) | imdilate(egl_strict, strel('disk', 40*scaling)); 
ml_layer_mask = ml_layer_mask > 0;  
ml_layer_mask = imfill(ml_layer_mask, 'holes'); 
deep_nuclei = all_cerebellum .* (1- ml_layer_mask); 
deep_nuclei = deep_nuclei > 0; 
% compensate for the initial dialtion 
deep_nuclei = imdilate(deep_nuclei, strel('disk', 40*scaling));
deep_nuclei = bwareafilt(deep_nuclei > 0,[10000 Inf]);
deep_nuclei_near_egl = im_get_overlap(deep_nuclei, egl_mask , 0.8, 100 );
deep_nuclei_near_egl = deep_nuclei_near_egl > 5; 
deep_nuclei = deep_nuclei .* (1 - deep_nuclei_near_egl); 

ml_layer_mask = (ml_layer_mask - igl_in - egl_mask - deep_nuclei) .* all_cerebellum;
ml_layer_mask = ml_layer_mask > 0;
ml_layer_mask = imopen(ml_layer_mask, strel('disk', 3*scaling)); 
ml_layer_mask = bwareafilt(ml_layer_mask > 0,[20000 Inf]);
pc_layer_mask = pc_layer_mask .* ml_layer_mask;

ml_layer_mask =  ml_layer_mask - pc_layer_bin;
ml_layer_mask = ml_layer_mask > 0;

dwl_final = dwl_final | deep_nuclei;
dwl_final = imdilate(dwl_final, strel('disk',40*scaling));
dwl_final = dwl_final .* (1-c_final);
dwl_final = imfill(dwl_final, 'holes');



igl_in =  (c_final | dwl_final ) .* all_cerebellum; 
igl_in = imfill(igl_in, 'holes');
dwl_final = igl_in - c_final; 
dwl_final = bwareafilt(dwl_final > 0,[200000 Inf]);
% refine the ML layer
ml_layer_mask = ml_layer_mask .* (1 - dwl_final);


ml_in = igl_in | ml_layer_mask;
% dilate to fill up to the iegl 
ml_in = imclose(ml_in, strel('disk',10*scaling));
ml_in = imfill(ml_in, 'holes');
ml_in = imdilate(ml_in, strel('disk',10*scaling));
ml_in = bwareafilt(ml_in > 0,[5000 Inf]);
ml_in = imdilate(ml_in, strel('disk',10*scaling));
%ml_in = imclose(ml_in, strel('disk',2));
ml_in = imfill(ml_in, 'holes');
ml_in = bwareafilt(ml_in > 0,[20000 Inf]);

ml_in = ml_in .* (1- a_bin); 
ml_in = ml_in == 1; 
ml_in = bwareafilt(ml_in > 0,[200000 Inf]);
ml_in = imfill(ml_in, 'holes');


% handle the case that a piece of the egl is disrupted and connects the ML
% to the CSF 
%imshow(imoverlay(egl_mask_large, ml_in, [0 1 0]))
egl_mask_large = imdilate(egl_mask, strel('disk', 20*scaling));
egl_mask_large = imfill(egl_mask_large, 'holes'); 

b_ml_intersect = im_get_overlap(egl_mask_large, ml_in , 0.8, 100 );
b_ml_intersect = b_ml_intersect >= 9; 

% release the pieces that were pre
%b_ml_intersect
%intersect_labels = zeros(size(all_cerebellum));
%intersect_labels(egl_mask == 1) = 1;
%intersect_labels(a_bin == 1) = 2;
%intersect_labels(ml_in == 1) = 3;
%imagesc(intersect_labels)
%(all_cerebellum .* 1) + (a_bin .* 2)  + (ml_in .* 3)  + (b_ml_intersect .* 4); 

ml_in = ml_in .* (1-b_ml_intersect); 
ml_in = ml_in == 1;
imshow(imoverlay(a_bin | igl_in, bwperim(ml_in), [0 1 0]))

%redilate to the iegl 
ml_in = imclose(ml_in, strel('disk',5*scaling));
ml_in = imfill(ml_in, 'holes');
ml_in = imdilate(ml_in, strel('disk',5*scaling));
ml_in = bwareafilt(ml_in > 0,[50000 Inf]);
ml_in = ml_in - a_bin - b_bin; 
ml_in = ml_in == 1;
ml_in = imfill(ml_in, 'holes');


%resubtract the oegl


% segregate the oegl 
oegl_in = ml_in | a_bin; 
%imshow(imoverlay(a_bin, b_bin, [0 1 0]))
%imshow(imoverlay(oegl_out, ml_bin, [0 1 0]))

ml_in = all_cerebellum - oegl_out;


% need to preserve the csf and pial boundary 
iegl_out = (1 - oegl_in) .* (all_cerebellum);
iegl_out = iegl_out  > 0 ;
iegl_out = bwareafilt(iegl_out > 0,[50 Inf]);
iegl_out = iegl_out | b_bin;
iegl_out = bwareafilt(iegl_out > 0,[100 Inf]);
iegl_out =  iegl_out - oegl_in;
iegl_out = iegl_out  > 0 ;
iegl_out = imclose(iegl_out, strel('disk',10*scaling));
iegl_out = imdilate(iegl_out, strel('disk',10*scaling));
iegl_out  = iegl_out - a_bin; 
iegl_out = iegl_out  > 0 ;
iegl_out_mask = imdilate(iegl_out, strel('disk',20*scaling));

% now dilate the oegl  
oegl_in = bwareafilt(oegl_in > 0,[500 Inf]);
oegl_in  = smallfill( oegl_in, 3000 );
oegl_in = imclose(oegl_in, strel('disk',10*scaling));
oegl_in  = smallfill( oegl_in, 500 );
oegl_in_mask = imdilate(oegl_in, strel('disk',20*scaling));

% now intersect the dilated images 
egl_intersect = im_get_overlap(oegl_in, iegl_out , 0.5, 100 );
egl_intersect = egl_intersect  > 4; 

oegl_in = ml_in | a_bin| egl_intersect;
oegl_in = bwareafilt(oegl_in > 0,[200 Inf]);
oegl_in  = smallfill( oegl_in, 3000 );
iegl_out = iegl_out - egl_intersect; 
iegl_out = iegl_out == 1;

% very sligh dilation 
oegl_in = imdilate(oegl_in, strel('disk',2));

a_final = a_bin .* oegl_in .* (1- ml_in) .* (1-igl_in) ; % strips the temporary G0 arrested cells in iEGL
b_final = b_bin  .* (1- ml_in) .* (1-igl_in) ;  %.* (1- oegl_in);
b_final = b_final == 1;
%b_final = imdilate(b_final, strel('disk',2));

ml_final = ml_in - igl_in;

pc_bin_filt = pc_bin .* ml_final;

%imshow(imoverlay((igl_in | a_bin | ml_layer_mask), b_bin, [0 1 0]))
%egl_col(:,:,1) = double(bwperim(a_bin .* egl_mask));
%egl_col(:,:,2) = double(bwperim(b_bin .* egl_mask));
%egl_col(:,:,3) = a_nfg_h;
%imshow(uint8(255*mat2gray(egl_col)));
    
clear  b_nfg_l c_nfg_l cerebellum_boundary boundary_y boundary_x

figure

montage([a_nf b_nf c_nf; ...
    a_final b_final c_final;...
    ml_final dwl_final all_outline; ...
    pc_bin_filt pc_layer_bin all_cerebellum])


MyMontage = getframe(gca);
imwrite(MyMontage.cdata,[file(1:(end-4)) '_montage.tif'],'tif');

%% Plot the purkinje cells  

%pc_col(:,:,1) = pc_bin;
%pc_col(:,:,2) = pc_mask_out;
%pc_col(:,:,3) = pc_mask_in;
%imshow(uint8(255*mat2gray(pc_col)));
%imshow(imoverlay(a, bwperim(  pc_bin_filt)))


%% Refine norm masks by setting non debatable regions to 1 and 0 for others 


set_bin = zeros(size(all_cerebellum));
set_bin(all_cerebellum == 1) = 1;

set_bin(ml_final == 1) = 5;
set_bin(pc_layer_bin == 1) = 7;
set_bin(dwl_final == 1) = 6;
set_bin(c_final == 1) = 4;
set_bin(b_final == 1) = 3;
set_bin(a_final == 1) = 2;
set_bin(pc_bin_filt == 1) = 7;
%set_bin(pia_outline == 1) = 0;
figure
imagesc(set_bin)

col_map =  [1 1 1
            0.9333333 0.1725490 0.1725490
            0.0000000 0.8039216 0.4000000                        
            0.9333333 0.5098039 0.3843137            
            1 1 1
            1 1 1 
            0.3607843 0.6745098 0.9333333];
show_bin = set_bin; 
b_inc = imdilate(b_final, strel('disk', 8*scaling));
show_bin(b_inc == 1) = 3;      
show_bin(ml_final == 1) = 5;
show_bin(a_final == 1) = 2;
show_bin(pc_bin_filt == 1) = 7;
show_bin(pc_layer_bin == 1) = 7;

imshow(label2rgb(show_bin, col_map, 'w'))            



imwrite(label2rgb(set_bin, 'jet', 'k'),[file(1:(end-4)) '_segments.tif'])

imwrite(double(a) .* double(a_bin),[file(1:(end-4)) '_a_iegl_only.tif'])
imwrite(double(a_nf) .* double(a_bin),[file(1:(end-4)) '_anf_iegl_only.tif'])
imwrite(imoverlay(a_nf, bwperim(a_final), [1 0 0]),[file(1:(end-4)) '_anf_overlay.tif'])
imwrite(imoverlay(mat2gray(a), bwperim(a_final), [1 0 0]),[file(1:(end-4)) '_a_overlay.tif'])
imwrite(imoverlay(b_nf, bwperim(a_final), [1 0 0]),[file(1:(end-4)) '_bnf_overlay.tif'])
imwrite(imoverlay(mat2gray(b), bwperim(a_final), [1 0 0]),[file(1:(end-4)) '_b_overlay.tif'])


% get the threshold that has 50% markig of IGL and then measure the other
% layers
%a_set_thresh =  a_set_thresh > get_thresh_by_area( a_set_thresh, c_final, 0.5 );
%b_set_thresh = b_set_thresh > get_thresh_by_area( b_set_thresh, c_final, 0.5 );


%imshow(imoverlay(a_nf,  c_final, [0 1 0]))
N = max(max(set_bin));

layer_vals = regionprops(set_bin,'PixelIdxList','PixelList', 'centroid', 'Area');


image_cell = {a_nf, b_nf, c_nf,a, b,  a_bin_1, b_bin_25}; 

sum_im_cell = repmat({zeros(1,N)},1,size(image_cell,2));
%median_im_cell = repmat({zeros(1,N)},1,size(image_cell,2));
mean_im_cell = repmat({zeros(1,N)},1,size(image_cell,2));



% need to vectorize this nested loop 
for n=1:length(layer_vals)  %for each cell in the image
    for i=1:size(image_cell, 2)
        sum_im_cell{i}(n) = sum(image_cell{i}(layer_vals(n).PixelIdxList));
        %median_im_cell{i}(n) = median(image_cell{i}(layer_vals(n).PixelIdxList));
        mean_im_cell{i}(n) = mean(image_cell{i}(layer_vals(n).PixelIdxList));
    end 
end

T = table([layer_vals.Area]', [mean_im_cell{1}]', [mean_im_cell{2}]', [mean_im_cell{3}]', [mean_im_cell{4}]', [mean_im_cell{5}]', [mean_im_cell{6}]', [mean_im_cell{7}]');
T.Properties.VariableNames = {'pixel_area', 'p27_mean', 'DAPI_mean', 'NeuN_mean', 'p27_unfilt', 'DAPI_unfilt', 'p27_thresh_mean', 'DAPI_thresh_mean' };

% 1 um = 0.9767 pixels 
% 1 mm = 976.7 pixels 
% 1 mm^2 = 9.5394e+05 pixels
mms2pix = 976.7^2;
T.area =  T.pixel_area / mms2pix;
region_names = {'all_cerebellum', 'iEGL', 'oEGL', 'IGL', 'ML', 'DWL', 'PCL'};  
T.region_names = region_names';

T1 = T;


set_bin = set_bin .* (1 - all_outline_mask);

layer_vals = regionprops(set_bin,'PixelIdxList','PixelList', 'centroid', 'Area');


image_cell = {a_nf, b_nf, c_nf,a, b,  a_bin_1, b_bin_25}; 

sum_im_cell = repmat({zeros(1,N)},1,size(image_cell,2));
%median_im_cell = repmat({zeros(1,N)},1,size(image_cell,2));
mean_im_cell = repmat({zeros(1,N)},1,size(image_cell,2));



% need to vectorize this nested loop 
for n=1:length(layer_vals)  %for each cell in the image
    for i=1:size(image_cell, 2)
        sum_im_cell{i}(n) = sum(image_cell{i}(layer_vals(n).PixelIdxList));
        %median_im_cell{i}(n) = median(image_cell{i}(layer_vals(n).PixelIdxList));
        mean_im_cell{i}(n) = mean(image_cell{i}(layer_vals(n).PixelIdxList));
    end 
end

T = table([layer_vals.Area]', [mean_im_cell{1}]', [mean_im_cell{2}]', [mean_im_cell{3}]', [mean_im_cell{4}]', [mean_im_cell{5}]', [mean_im_cell{6}]', [mean_im_cell{7}]');
T.Properties.VariableNames = {'pixel_area', 'p27_mean', 'DAPI_mean', 'NeuN_mean', 'p27_unfilt', 'DAPI_unfilt', 'p27_thresh_mean', 'DAPI_thresh_mean' };

% 1 um = 0.9767 pixels 
% 1 mm = 976.7 pixels 
% 1 mm^2 = 9.5394e+05 pixels
mms2pix = 976.7^2;
T.area =  T.pixel_area / mms2pix;
region_names = {'all_cerebellum', 'iEGL', 'oEGL', 'IGL', 'ML', 'DWL', 'PCL'};  
T.region_names = region_names';

T2 = T;

T1.inner = repmat({'full'}, size(T1,1),1);
T2.inner = repmat({'inner'}, size(T2,1),1);

T = [T1; T2];

layer = T;

return
end

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


