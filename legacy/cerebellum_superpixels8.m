function layer = cerebellum_superpixels8(file)  

a = imread(file,1);
b = imread(file, 2);
c = imread(file, 3);

%a = illumination_cor(a);
%b = illumination_cor(b);
%c = illumination_cor(c);
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
[a_ps_nmxs,location] = nonmaxsup(a_phasesym, orientation, 3);
disp('phase symmetry done')

col(:,:,1) = a_nf;
col(:,:,3) = b_nf;
col(:,:,2) = c_nf;

% remove low intensity parts of image 
a_nfg_h = background_subtract(mat2gray(a_nfg, [0.4 1])); 
b_nfg_h = background_subtract(mat2gray(b_nfg, [0.4 1]));
c_nfg_h = background_subtract(mat2gray(c_nfg, [0.4 1]));
a_nfg_l = mat2gray(a_nfg, [0 0.4]); 
b_nfg_l = mat2gray(b_nfg, [0 0.4]);
c_nfg_l = mat2gray(c_nfg, [0 0.4]);

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
pc_bin = imopen(pc_bin, strel('disk',2));
pc_bin = bwareafilt(pc_bin > 0,[40 Inf]);
pc_bin = imdilate(pc_bin,strel('disk', 4));

%imshow(label2rgb(keeperBlobsImage, 'winter', 'k'))
%imshow(imoverlay(a, bwperim(  pc_bin)))

clear measurements allAreas allPerimeters circularities roundObjectsIndexes
clear keeperBlobsImage
%clear a_nf b_nf c_nf a b c 




all_cerebellum = b_nfad  + a_nfad;
all_cerebellum = background_subtract(all_cerebellum);
all_cerebellum = all_cerebellum > 0.1; 
% handles the case when there is a breach in the EGL 
all_cerebellum_large = imclose(all_cerebellum, strel('disk', 20));
all_cerebellum_large = imfill(all_cerebellum_large, 'holes');
all_cerebellum_large = all_cerebellum_large - all_cerebellum;
all_cerebellum_large = imclose(all_cerebellum_large, strel('disk', 5));
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
all_cerebellum = imclose(all_outline, strel('disk', 10));
all_cerebellum = imfill(all_cerebellum, 'holes');
all_outline = bwperim(all_cerebellum);
all_outline_mask = imdilate(all_outline,strel('disk', 50));
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
c_bin = thresh_by_area( c_norm, all_cerebellum, 0.3 );
c_bin = c_bin .* (1 - all_outline_mask);
c_bin = bwareafilt(c_bin > 0,[5000 Inf]);
c_bin = imclose(c_bin, strel('disk',5));
c_mask = imclose(c_bin, strel('disk',50));
c_mask = imdilate(c_mask, strel('disk',8));

c_whole = im_get_overlap(c_mask, c_nfg_h > 0.4, 0.1, 200);
ca_whole = im_get_overlap(c_mask, a_nfg_h > 0.5, 0.6, 200);
cb_whole = im_get_overlap(c_mask, b_nfg_h > 0.4, 0.1, 200);
c_final = (ca_whole + c_whole + cb_whole)  > 0 ; 
c_final = imclose(c_final, strel('disk', 10));
c_final = imfill(c_final, 'holes');


disp('layer norm done')


% get the white matter areas deep to the igl 
% they have relatively high levels of p27 comapred to the ML 
dwl_norm = (double(a_nfg_l) - double(c_nfg_l) - double(a_nfg_h) - (double(a_nfg_h)))  ; 
dwl_norm = mat2gray(dwl_norm) .* all_cerebellum;
dwl_norm = ordfilt2(dwl_norm,25,ones(5,5));
dwl_norm =  imgaussfilt(dwl_norm,5);
dwl_norm = background_subtract(dwl_norm);


dwl_bin = thresh_by_area( dwl_norm, all_cerebellum, 0.08 );
dwl_bin = dwl_bin .* (1 - imdilate(all_outline_mask,strel('disk', 70)));
dwl_bin = dwl_bin .* (1 - c_bin);

dwl_mask = imdilate(dwl_bin, strel('disk',2));
dwl_mask = imclose(dwl_mask, strel('disk',15));
dwl_mask = imdilate(dwl_mask, strel('disk',10));
dwl_mask = bwareafilt(dwl_mask > 0,[4000 Inf]);

igl_in = c_bin | c_final | dwl_bin ;
dwl_phase = (a_phasesym .* igl_in > 0.4);
igl_in = igl_in | dwl_phase ;

igl_in = imdilate(igl_in, strel('disk',10));
igl_in = imclose(igl_in, strel('disk',10));
igl_in = imfill(igl_in, 'holes');
igl_in = bwareafilt(igl_in > 0,[70000 Inf]);


% TO DO, needs to be more robust to pieces broken off DWL
dwl_final = igl_in - c_final;
dwl_final = dwl_final | dwl_phase;
dwl_final = imerode(dwl_final, strel('disk',10));
dwl_final = imdilate(dwl_final, strel('disk',20));
dwl_final = bwareafilt(dwl_final > 0,[500000 Inf]);
igl_in = c_final | dwl_final; 
dwl_final = igl_in - c_final; 


% Robustly call the EGL 
%egl_norm = (double(b_nfad) + double(a_nfad)) .* (1- igl_in) .* all_cerebellum;
egl_norm = double(b_nfg)  .* (1- igl_in) .* all_cerebellum;
egl_norm = ordfilt2(egl_norm,25,ones(5,5));
egl_norm = background_subtract(mat2gray(egl_norm));
egl_bin = thresh_by_area( egl_norm, all_cerebellum, 0.15 );
egl_mask = imclose(egl_bin, strel('disk',5));
egl_mask = imclose(egl_mask, strel('disk',5));
egl_mask = bwareafilt(egl_mask > 0,[5000 Inf]);
egl_phase = im_get_overlap(egl_mask, a_phasesym > 0.1, 0.1, 200 );
egl_mask = imclose(egl_mask | egl_phase, strel('disk',10));
egl_mask = imdilate(egl_mask, strel('disk',10));
egl_mask = imerode(egl_mask, strel('disk',5));
egl_mask = imdilate(egl_mask, strel('disk',10));
egl_mask = imclose(egl_mask | egl_phase, strel('disk',10));


%egl_thin = bwmorph(egl_mask,'thin', Inf);
%egl_thin_mask = imdilate(egl_thin, strel('disk',10));

% refine the purkinje cell zone, which includes more then PC 
pc_layer_mask = imdilate(igl_in, strel('disk', 40));
pc_layer_mask = pc_layer_mask - igl_in - egl_mask;
pc_layer_mask = pc_layer_mask == 1;
pc_layer_bin = (pc_layer_mask .* b_nfg) > 0.4;

ml_layer_mask = imdilate(igl_in, strel('disk', 50)); 
ml_layer_mask = ml_layer_mask + imdilate(egl_mask, strel('disk', 50)); 
ml_layer_mask = imfill(ml_layer_mask, 'holes'); 
deep_nuclei = all_cerebellum - ml_layer_mask; 
% compensate for the initial dialtion 
deep_nuclei = imdilate(deep_nuclei, strel('disk', 40));
ml_layer_mask = (ml_layer_mask - igl_in - egl_mask - deep_nuclei) .* all_cerebellum;
ml_layer_mask = ml_layer_mask > 0;
ml_layer_mask = imopen(ml_layer_mask, strel('disk', 3)); 
ml_layer_mask = bwareafilt(ml_layer_mask > 0,[20000 Inf]);
pc_layer_mask = pc_layer_mask .* ml_layer_mask;

ml_layer_mask =  ml_layer_mask - pc_layer_bin;
ml_layer_mask = ml_layer_mask > 0;

dwl_final = dwl_final | deep_nuclei;
dwl_final = imdilate(dwl_final, strel('disk',20));
dwl_final = bwareafilt(dwl_final > 0,[1000000 Inf]);
dwl_final = imfill(dwl_final, 'holes');
igl_in =  (c_final | dwl_final ) .* all_cerebellum; 
dwl_final = igl_in - c_final; 

% refine the ML layer
ml_layer_mask = ml_layer_mask .* (1 - dwl_final);



a_nfad_mean = mean(mean(a_nfad));
a_nfad_egl = a_nfad;
a_nfad_egl((1-all_cerebellum) > 0) = a_nfad_mean;
T = adaptthresh(a_nfad_egl , 0.35,'NeighborhoodSize' , 2*floor(size(a_nfg)/25)+1);
a_bin = imbinarize(a_nfad_egl ,T);
a_bin = a_bin .* egl_mask;
%imshow(imoverlay(a_nf, bwperim(a_bin ),  [0 1 0]))

b_nfad_mean = mean(mean(b_nfad));
b_nfad_egl = b_nfad;
b_nfad_egl((1-all_cerebellum) > 0) = b_nfad_mean;
T = adaptthresh(b_nfad_egl , 0.6,'NeighborhoodSize' , 2*floor(size(b_nfad_egl)/25)+1);
b_bin = imbinarize(b_nfad_egl ,T);
b_bin = b_bin .* egl_mask;

%imshow(imoverlay(b_nf, bwperim(b_bin ),  [0 1 0]))
%imshow(imoverlay(b_nf, bwperim(egl_mask ),  [0 1 0]))

b_bin = b_bin - a_bin; 
b_mask = imclose(b_bin, strel('disk', 10));
b_mask = imdilate(b_bin, strel('disk', 20));
b_mask = b_mask - a_bin; 
b_mask = b_mask == 1;



ml_in = igl_in | ml_layer_mask;
% dilate to fill up to the iegl 
ml_in = imclose(ml_in, strel('disk',10));
ml_in = imfill(ml_in, 'holes');
ml_in = imdilate(ml_in, strel('disk',10));
ml_in = bwareafilt(ml_in > 0,[5000 Inf]);
ml_in = imdilate(ml_in, strel('disk',10));
%ml_in = imclose(ml_in, strel('disk',2));
ml_in = imfill(ml_in, 'holes');
ml_in = bwareafilt(ml_in > 0,[20000 Inf]);

ml_in = ml_in - a_bin - b_bin; 
ml_in = ml_in == 1; 
ml_in = bwareafilt(ml_in > 0,[100000 Inf]);


% handle the case that a piece of the egl is disrupted and connects the ML
% to the CSF 
%imshow(imoverlay(egl_mask_large, ml_in, [0 1 0]))
egl_mask_large = imdilate(egl_mask, strel('disk', 20));
egl_mask_large = imfill(egl_mask_large, 'holes'); 

b_ml_intersect = im_get_overlap(egl_mask_large, ml_in , 0.8, 100 );
b_ml_intersect = b_ml_intersect >= 9; 
ml_in = ml_in - b_ml_intersect; 
ml_in = ml_in == 1;
imshow(imoverlay(a_bin, ml_in, [0 1 0]))

%redilate to the iegl 
ml_in = imclose(ml_in, strel('disk',5));
ml_in = imfill(ml_in, 'holes');
ml_in = imdilate(ml_in, strel('disk',5));
ml_in = bwareafilt(ml_in > 0,[50000 Inf]);
ml_in = ml_in - a_bin - b_bin; 
ml_in = ml_in == 1;
%resubtract the oegl

% segregate the oegl 
oegl_in = ml_in | a_bin; 
% need to preserve the csf and pial boundary 
iegl_out = (1 - oegl_in) .* (all_cerebellum);
iegl_out = iegl_out  > 0 ;
iegl_out = bwareafilt(iegl_out > 0,[50 Inf]);
iegl_out = iegl_out | b_bin;
iegl_out = bwareafilt(iegl_out > 0,[100 Inf]);
iegl_out =  iegl_out - oegl_in;
iegl_out = iegl_out  > 0 ;
iegl_out = imclose(iegl_out, strel('disk',10));
iegl_out = imdilate(iegl_out, strel('disk',10));
iegl_out  = iegl_out - a_bin; 
iegl_out = iegl_out  > 0 ;

% now dilate the oegl  
oegl_in = imclose(oegl_in, strel('disk',10));

% now intersect the dilated images 
egl_intersect = im_get_overlap(oegl_in, iegl_out , 0.5, 100 );
egl_intersect = egl_intersect  > 4; 

oegl_in = ml_in | a_bin| egl_intersect;
iegl_out = iegl_out - egl_intersect; 
iegl_out = iegl_out == 1;

% very sligh dilation 
oegl_in = imdilate(oegl_in, strel('disk',2));

a_final = a_bin;
b_final = b_bin - oegl_in;

b_final = b_final == 1;
b_final = imdilate(b_final, strel('disk',2));
ml_final = oegl_in -(a_bin | c_final  | dwl_final | b_final | pc_layer_bin );



imshow(imoverlay((igl_in | a_bin | ml_layer_mask), b_bin, [0 1 0]))

egl_col(:,:,1) = double(bwperim(a_bin .* egl_mask));
egl_col(:,:,2) = double(bwperim(b_bin .* egl_mask));
egl_col(:,:,3) = a_nfg_h;
imshow(uint8(255*mat2gray(egl_col)));
    
clear a_nfg_l b_nfg_l c_nfg_l cerebellum_boundary boundary_y boundary_x



%montage([a_nf b_nf c_nf ...
%    ml_norm dwl_norm all_outline;...
%    a_bin b_bin c_bin; ...
%    ml_bin dwl_bin all_cerebellum])

%% -----------------------------------------------------------------------
% FIND THE PIAL BOUNDARY 
% identifying the pial boundary is a little tricky 
% the position has to exact in folia, dividing the outerEGL using the
% skeleton isn't perfect
% using phase symmetry we can enrich the line enough to get a signal but if
% the layers aren't exactly opposed the signal drops 
% use a combination of approaches

% use EGL channels to define the approximate location of the pia within the
% buried portion of the folia 
%pia = imbinarize(mat2gray(b_norm+a_norm-c_norm-ml_norm-dwl_norm), 'adaptive', 'Sensitivity',0.8);
%pia_norm =egl_mask;
%pia_norm(pia_norm < 0) = 0; 
%pia_bin = thresh_by_area( pia_norm, all_cerebellum, 0.2 );

pia = imdilate(egl_mask, strel('disk',15));
pia = imerode(pia, strel('disk',10));;
pia = bwareafilt(pia > 0,[10000 Inf]);
pia = imdilate(pia, strel('disk',10));
pia = imerode(pia, strel('disk',10));
pia_thin = bwmorph(pia,'thin', Inf);




% refine the iegl or oegl mask 
iegl_refine = (1-(all_cerebellum - a_all_mask ) ).* (pia);
oegl_refine = pia - iegl_refine ; 
oegl_refine = imerode(oegl_refine, strel('disk',1));
oegl_refine = bwareafilt(oegl_refine > 0,[5000 Inf]);
oegl_deep = oegl_refine .*  (1- all_outline_mask) .* pia;

iegl_in = iegl_in .* (1- oegl_deep);



se = strel('disk',15);
pia_mask = imdilate(pia_thin, se);
% -------------------------------------------------------------------------
% Use phase symmetry to deliniate the pial border 
pia_phasesym = ((b_phasesym + a_phasesym) .* pia_mask);
[pia_line,location] = nonmaxsup(pia_phasesym, orientation, 3);
% TO DO, need to tune this parameter 
pia_line = pia_line > 0.5; 
pia_line = bwareafilt(pia_line > 0,[10 Inf]);
%pia_bin2 = pia_phasesym >0.5 ;

% ------------------------------------------------------------------------
% Use the csf space near the pia to make a more clear line 
csf_space = ((1- b_nfad)  > 0.6) .* all_cerebellum .* (1 - (a_bin));
csf_space = csf_space .* pia_mask ;

% -----------------------------------------------------------------------
% Try and define the a mask that has the inner EGL and contained contents
iegl_in = igl_in | ml_layer_mask | a_bin;

% subtract the csf from the egl_in mask 
iegl_in = iegl_in  - csf_space;

iegl_in = imclose(iegl_in, strel('disk',3));
iegl_in = bwareafilt(iegl_in > 0 ,[10000 Inf]);



% define the center of the deep egl segments
pia = pia .* (1-a_all_mask); 
pia_center_thin = bwmorph(pia,'thin', Inf);
pia_center_thin = pia_center_thin .* (1- iegl_in);
pia_center_thin = pia_center_thin .* (1- all_outline_mask);

% add the csf to the pia mask  
csf_edge = edge(csf_space, 'canny');
csf_edge = imdilate(csf_edge, strel('disk',3));
csf_edge = imerode(csf_edge, strel('disk',2));
csf_edge = imfill(csf_edge,'holes');
csf_edge = bwmorph(csf_edge,'skel', Inf);


pia_inner = pia_line; %.* (1-all_outline_mask);
% strip segments overlaping with outer EGL
%pia_inner_bwl = bwlabel(pia_inner);
%pia_inner = ~ismember(pia_inner_bwl , unique((a_all_mask | c_mask ) .* pia_inner_bwl));
pia_inner = imdilate(pia_inner, strel('disk',3));
csf_edge = imdilate(csf_edge, strel('disk',3));
pia_center_thin = imdilate(pia_center_thin, strel('disk',3));
pia_sum = (pia_inner + csf_edge + pia_center_thin) .* (1-all_outline_mask); 
[pia_sum_nms,location] = nonmaxsup(pia_sum, orientation, 3);
pia_sum = pia_sum_nms > 1; 



%pia_col(:,:,1) = pia_inner;
%pia_col(:,:,2) = csf_edge;
%pia_col(:,:,3) = pia_center_thin;
%imshow(uint8(255*mat2gray(pia_col))); 

pia_inner = bwareafilt(pia_sum > 0,[3 Inf]);

%pia_inner = bwmorph(pia_inner_mask,'thin', Inf);
%pia_inner = pia_inner .* (1 - iegl_in);
pia_inner_fill = filledgegaps(pia_inner, 3);
pia_inner_fill = filledgegaps(pia_inner_fill, 5);
pia_inner_fill = pia_inner_fill .* (1 - (ml_in)) .* (1-all_outline_mask);
pia_inner_fill = bwareafilt(pia_inner_fill > 0,[10 Inf]);
pia_inner_fill = filledgegaps(pia_inner_fill, 11);
pia_inner_fill = filledgegaps(pia_inner_fill, 21);
pia_inner_fill = pia_inner_fill .* (1 - (ml_in)) .* (1-all_outline_mask);
pia_inner_fill = bwareafilt(pia_inner_fill > 0,[20 Inf]);
%pia_inner_fill = bwareafilt(pia_inner_fill > 0,[20 Inf]);
pia_inner_fill = filledgegaps(pia_inner_fill, 41);
%pia_inner_fill = pia_inner_fill .* (1 - (a_norm > 0.2)) .* (1-all_outline_mask);
pia_inner_fill = filledgegaps(pia_inner_fill, 61);
%pia_inner_fill = pia_inner_fill .* (1 - (a_norm > 0.2)) .* (1-all_outline_mask);
pia_inner_fill = filledgegaps(pia_inner_fill, 81);
%pia_inner_fill = pia_inner_fill .* (1 - (a_norm > 0.2)) .* (1-all_outline_mask);
pia_inner_fill = filledgegaps(pia_inner_fill, 101);
%montage([imoverlay(a_norm, pia_inner, [0 1 0]) imoverlay(a_norm, pia_inner_fill, [0 1 0])])
%pia_inner_fill = pia_inner_fill .* (1 - (a_norm > 0.2)) .* (1-all_outline_mask);
pia_inner_fill = filledgegaps(pia_inner_fill, 121);
pia_inner_fill = filledgegaps(pia_inner_fill, 141);
pia_inner_fill = filledgegaps(pia_inner_fill, 161);
%pia_inner_fill = pia_inner_fill .* (1 - (a_norm > 0.2)) .* (1-all_outline_mask);
pia_inner_fill = pia_inner_fill .* (1 - (ml_in)) .* (1-all_outline_mask);
pia_inner = pia_inner_fill;
pia_deep = pia_inner .* (1- ml_in); 

pia_deep = bwareafilt(pia_deep > 0,[100 Inf]);
%imshow(imoverlay(b_nfad, pia_inner_fill, [0 1 0 ] ) )
%imshow(imoverlay(pia_sum_nms, pia_inner_fill, [0 1 0 ] ) )

imshow(pia_deep | pia_sum)

%all_outline .* 1 - (a_norm > 0.2)


% dilate the pia as marked by phase sym 
pia_inner_mask = imdilate(pia_deep, strel('disk',10));
pia_inner_mask = pia_inner_mask .* all_cerebellum;


% calculate overlap with the csf space 
csf_lab = bwlabel(csf_space);
%imshow(label2rgb(csf_lab))
all_pia_ov = unique(csf_lab .* pia_inner_mask);
all_pia_ov = sort(all_pia_ov);
all_pia_ov = all_pia_ov(2:end);
pia_ov_t = histc(reshape(csf_lab .* pia_inner_mask,[],1), all_pia_ov) ./ histc(reshape(csf_lab ,[],1), all_pia_ov);

% keep high overlaping CSF space 
pia_ov_sel = all_pia_ov(pia_ov_t > 0.001);
pia_csf_ov = ismember(csf_lab, pia_ov_sel );

% for every pia_deep_segment get the start and end and find the min
% distance to the the pia outline then burn a line to this point 
pia_deep_bwl = bwlabel(pia_deep);  
pia_discard = unique(pia_deep_bwl .* c_mask);
junctions = join_masks(pia_deep_bwl, all_outline, pia_discard);
pia_deep(ismember(pia_deep_bwl, pia_discard)) = 0;
pia_deep = pia_deep | junctions;
pia_outline = all_outline | pia_deep | ismember(csf_lab, pia_ov_sel ) ;

clear csf_edge pia_line pia_phasesym
clear pia_deep_bwl pia_discard junctions pia_csf_ov pia_ov_sel pia_ov_t
clear all_pia_ov csf_lab pia_inner_mask pia_inner_fill pia_sum_nms pia_sum




%% refine the purkinje cells  


pc_mask_out = imerode(ml_in, strel('disk', 50));
pc_mask_in = imerode(ml_in, strel('disk', 160));

pc_col(:,:,1) = pc_bin;
pc_col(:,:,2) = pc_mask_out;
pc_col(:,:,3) = pc_mask_in;
%imshow(uint8(255*mat2gray(pc_col)));
%

pc_bin_filt = pc_bin .* pc_mask_out .* (1-pc_mask_in);

%imshow(imoverlay(a, bwperim(  pc_bin_filt)))


%% Refine norm masks by setting non debatable regions to 1 and 0 for others 


set_bin = zeros(size(all_cerebellum));
set_bin(all_cerebellum == 1) = 1;
set_bin(ml_final == 1) = 5;
set_bin(c_final == 1) = 4;
set_bin(pc_layer_bin == 1) = 7;
set_bin(dwl_final == 1) = 6;
set_bin(b_final == 1) = 3;
set_bin(a_final == 1) = 2;

%set_bin(pc_bin_filt == 1) = 7;
%set_bin(pia_outline == 1) = 0;
imagesc(set_bin)

%imshow(imoverlay(a_nf,  c_final, [0 1 0]))


%% Make pseudocolor image for segmentation usnig superpixels 



cmyk_col(:,:,1) = a_norm_m .* (1 - (pia_outline | pc_bin_filt)) .* all_cerebellum;
cmyk_col(:,:,2) = b_norm_m .* (1 - (pia_outline )) .* all_cerebellum;
cmyk_col(:,:,3) = (dwl_norm_m + ml_norm_m) .* (1 - (pia_outline  | pc_bin_filt)) .* all_cerebellum;
cmyk_col(:,:,4) = c_norm_m .* (1 - (pia_outline )) .* all_cerebellum;

%inprof = iccread('USSheetfedCoated.icc');
%outprof = iccread('sRGB.icm');
%C = makecform('icc',inprof,outprof);
C = makecform('cmyk2srgb'); 
I_rgb = applycform(cmyk_col,C);

figure
imshow(I_rgb)




%% Perform superpixel segmentation 




[L,N] = superpixels(I_rgb,20000);
BW = boundarymask(L);
imshow(imoverlay(I_rgb,BW))

[L_col,N] = superpixels(col,20000);
BW = boundarymask(L_col);
imshow(imoverlay(col,BW))




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






