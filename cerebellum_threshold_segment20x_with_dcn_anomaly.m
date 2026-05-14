function layer = cerebellum_threshold_segment20x_with_dcn_anomaly(file, varargin)
% 20x with-DCN segmentation + anomaly detection (V1).
% Detects three classes of slide artifacts and excludes them from layer
% segmentation: bleaching (low all-channel signal), folds (high all-channel
% signal), stitching seams (linear low-signal strips). Detected anomalies
% become a new label class 9 in the output, and they are subtracted from
% all_cerebellum so that downstream layer thresholding doesn't see them.
%
% Outputs (added vs with_dcn version):
%   *_anomaly_bleach.tif  uint8 mask of bleaching detections
%   *_anomaly_fold.tif    uint8 mask of fold detections
%   *_anomaly_seam.tif    uint8 mask of seam detections
%   *_anomaly_all.tif     uint8 union mask
% set_bin label 9 = anomaly. col_map / region_names extended.

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


%% ANOMALY DETECTION V4: tear + bleach via combined-channel + flat-field
%
% Key insight (validated on s1_2 gold standard tear ground truth):
%   - Real tears have DAPI median ~0.10 (NOT zero) and combined-channel
%     median ~0.085. About half of normal tissue intensity.
%   - All 3 channels drop together (artifacts hit all 3 fluorophores;
%     biology has channel-specific signal — p27 in iEGL, NeuN in IGL/DCN).
%   - Tears have eccentricity > 0.95 and aspect_ratio > 4 (long thin slits).
%   - ML inter-cell low-DAPI network has eccentricity 0.5-0.9 (not as
%     elongated as a real tear). Eccentricity is the discriminator.
%   - Bleaching patches are roundish/rectangular: ecc < 0.85, extent > 0.5.
%
% Pipeline:
%   1. Use the existing all_cerebellum mask (already morphologically closed
%      and hole-filled, so it INCLUDES the tear interior as "tissue").
%   2. Tissue-aware flat-field correction on each channel (sigma 100*scaling
%      = ~50 um at 20x). Removes large-scale illumination inhomogeneity
%      that creates fake "dropouts" on slides like s2_5.
%   3. combined_ff = mean of the 3 corrected channels.
%   4. Erode tissue 50*scaling px (25 um) to exclude cerebellum-edge
%      artifacts (FFT bleed at the tissue boundary).
%   5. Per-pixel threshold combined_ff < 0.10 → dropout candidates.
%   6. Connected components → shape filter:
%        TEAR:   area >= 1000 AND (eccentricity >= 0.95 OR axR >= 4)
%        BLEACH: area >= 2000 AND eccentricity <  0.85 AND extent > 0.5
%        REJECT (most ML/inhomogeneity blobs): otherwise
%
% V1 fold + seam detectors were dropped (fold was firing on EGL).
% V2 used local-max DAPI which was too conservative (zero detections).
% V3 was V2 + flat-field; still too conservative.
% V4 ditches local-max in favor of combined-channel + shape filter and
% works on the actual signal levels observed in real tears.
% Cast to logical: existing all_cerebellum is double (zeros() initialization
% propagated through imclose/imfill). We need logical for mask indexing.
all_cerebellum_orig = logical(all_cerebellum);

% IMPORTANT: build a MORE PERMISSIVE tissue mask just for anomaly detection.
% The pipeline's all_cerebellum uses imclose(disk(20*scaling)=40px) which is
% too small to close BIG tears (the s1_2 gold tear is ~150 px wide). So
% the tear interior gets EXCLUDED from all_cerebellum and we never see it.
%
% Solution: build tissue_perm with imclose(disk(150 px) = 75 um close +
% fill_holes. This includes the tear interior. We use tissue_perm only for
% anomaly detection — it doesn't affect the layer pipeline downstream.
tissue_perm = imclose(all_cerebellum_orig, strel('disk', 150));
tissue_perm = imfill(tissue_perm, 'holes');

% Normalize each channel to [0,1] for thresholding
a_n = mat2gray(double(a));   % p27
b_n = mat2gray(double(b));   % DAPI (channel b at 20x = page 3 of TIFF)
c_n = mat2gray(double(c));   % NeuN

% Tissue-aware flat-field correction
% bg_local(x,y) = locally-averaged intensity inside tissue at scale sigma.
% Then ch_ff = ch / bg_local * mean(bg_local). This rescales each region
% to compensate for illumination inhomogeneity.
% Use tissue_perm (includes tears) so FF reference doesn't see the tear as
% "missing tissue" and treat tear neighbors as background.
ff_sigma = 100*scaling;   % ~50 um at 20x — bigger than layers, smaller than image
mask_d = double(tissue_perm);
mask_blurred = imgaussfilt(mask_d, ff_sigma);
mask_blurred = max(mask_blurred, 0.05);   % avoid divide-by-zero near edge

bg_a = imgaussfilt(a_n .* mask_d, ff_sigma) ./ mask_blurred;
bg_b = imgaussfilt(b_n .* mask_d, ff_sigma) ./ mask_blurred;
bg_c = imgaussfilt(c_n .* mask_d, ff_sigma) ./ mask_blurred;

bg_a_mean = mean(bg_a(tissue_perm));
bg_b_mean = mean(bg_b(tissue_perm));
bg_c_mean = mean(bg_c(tissue_perm));

a_ff = max(0, min(1, a_n ./ max(bg_a, 0.02) * bg_a_mean));
b_ff = max(0, min(1, b_n ./ max(bg_b, 0.02) * bg_b_mean));
c_ff = max(0, min(1, c_n ./ max(bg_c, 0.02) * bg_c_mean));

combined_ff = (a_ff + b_ff + c_ff) / 3;

% Erode the PERMISSIVE tissue mask to exclude border-artifact false positives.
% Reduced from 50*scaling=100px (~50 um) to 25*scaling=50px (~25 um) to
% expose more of the tear interior to detection.
% Trade-off: slight increase in edge false positives near cerebellum boundary.
tissue_inner = imerode(tissue_perm, strel('disk', 25*scaling));

% Per-pixel dropout: combined intensity below threshold AND inside inner tissue.
% Threshold 0.13 (was 0.10) — gold tears have combined median ~0.085 but q75
% extends to 0.20. Raising from 0.10 to 0.13 captures peripheral tear pixels
% (~50% of gold tear) without firing on normal tissue (q5 of normal = 0.142).
dropout = (combined_ff < 0.13) & tissue_inner;
dropout = imclose(dropout, strel('disk', 2*scaling));   % bridge 1-2 px noise gaps
% Note: do NOT bwareafilt yet — we want all components for shape analysis

% Connected components + per-component shape stats
cc = bwconncomp(dropout);
stats = regionprops(cc, 'Area', 'MajorAxisLength', 'MinorAxisLength', ...
                       'Eccentricity', 'Extent');

tear_mask = false(size(dropout));
bleach_mask = false(size(dropout));
n_tear = 0; n_bleach = 0;
for k = 1:length(stats)
    area = stats(k).Area;
    ecc = stats(k).Eccentricity;
    major = stats(k).MajorAxisLength;
    minor = max(stats(k).MinorAxisLength, 0.5);
    axR = major / minor;
    extent = stats(k).Extent;
    if area >= 1000 && (ecc >= 0.95 || axR >= 4)
        tear_mask(cc.PixelIdxList{k}) = true;
        n_tear = n_tear + 1;
    elseif area >= 2000 && ecc < 0.85 && extent > 0.5
        bleach_mask(cc.PixelIdxList{k}) = true;
        n_bleach = n_bleach + 1;
    end
end

anomaly_mask = tear_mask | bleach_mask;

cereb_px = max(sum(all_cerebellum_orig(:)), 1);
fprintf(['ANOMALY V4: tears=%d (%.2f%%) bleach=%d (%.2f%%) total=%.2f%% ' ...
         '(of cereb)\n'], ...
    n_tear, 100*sum(tear_mask(:))/cereb_px, ...
    n_bleach, 100*sum(bleach_mask(:))/cereb_px, ...
    100*sum(anomaly_mask(:))/cereb_px);

% Save masks for inspection
imwrite(uint8(255*tear_mask),    [file(1:(end-4)) '_anomaly_tear.tif']);
imwrite(uint8(255*bleach_mask),  [file(1:(end-4)) '_anomaly_bleach.tif']);
imwrite(uint8(255*anomaly_mask), [file(1:(end-4)) '_anomaly_all.tif']);

% Subtract anomaly from cerebellum mask — pipeline runs on cleaned tissue
all_cerebellum = all_cerebellum_orig & ~anomaly_mask;
disp('anomaly detection done')


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

% IGL FIX: original used (ca_whole | c_whole | cb_whole) — pure OR meant a
% region with high p27 OR high DAPI inside the smoothed-NeuN-top-40% area
% became IGL even with little raw NeuN. PC and basket cells (NeuN-, p27+/-,
% DAPI+) sitting just inside the IGL boundary got incorrectly absorbed into
% IGL. Symptom user reported: "non-NeuN cells gathering near the PCL".
%
% New rule: IGL membership requires raw NeuN intensity directly. Either
%   (a) confirmed by NeuN signal (c_whole), OR
%   (b) confirmed by p27 AND raw NeuN above c_level (ca_whole + neun_strict).
% Drop pure-DAPI-confirmation (cb_whole) entirely — that was the worst
% offender, since DAPI is high in PCs and basket cells with no NeuN.
neun_strict = c_nf > c_level;
c_final = c_whole > 0 | (ca_whole > 0 & neun_strict);
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

% PCL FIX: original used ANY DAPI-bright pixel in this ring as PCL.
% Problem: basket and stellate interneurons (8-10 um cells) also have bright
% DAPI nuclei and live in the molecular layer adjacent to PCs — they were
% being mislabeled as PCL. Real Purkinje cell SOMATA are 15-25 um diameter
% (~30-50 px at 20x = ~700-2000 px area). Filter pc_layer_bin by component
% area to keep only PC-sized blobs.
pc_layer_seed = (pc_layer_mask .* b_nfg) > 0.4;
pc_layer_seed = bwareafilt(pc_layer_seed > 0, [round(500*scaling^2) round(8000*scaling^2)]);
pc_layer_bin = pc_layer_seed;


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



%% Detect Deep Cerebellar Nuclei (DCN): NeuN+ clusters within DWL
% Clusters of large NeuN+ neurons in the deep white matter tract; biologically
% distinct from surrounding WM. Currently lumped into DWL — we split them out.
%
% IMPORTANT: DCN must come from DWL only — never from EGL/ML/IGL. We enforce
% this by eroding DWL before assigning DCN, so morphological closing/filling
% during seed cleanup can never push the DCN region into adjacent layers.
dcn_final = false(size(dwl_final));
if any(dwl_final(:))
    % Pre-compute the safe interior of DWL — DCN can only live here.
    % Erosion radius (10 px at 20x) is the safety margin against the close+fill
    % below spilling into ML/EGL territory at the DWL boundary.
    dwl_interior = imerode(dwl_final, strel('disk', 10*scaling));

    neun_in_dwl = c_nf(dwl_interior);
    if numel(neun_in_dwl) > 100
        % Otsu within the deep DWL gives a stable threshold even when there's
        % no DCN at all (otsu just splits noise from noise).
        thr_otsu = graythresh(neun_in_dwl);
        % Stricter NeuN floor: must be well above WM background. Without this,
        % otsu on slides with no DCN happily grabs noise.
        thr_mean = mean(neun_in_dwl) + 1.0 * std(double(neun_in_dwl));
        thr_floor = 0.20;   % NeuN must be >= 20% of dynamic range (anywhere)
        thr = max([thr_otsu, thr_mean, thr_floor]);

        % Seed: NeuN+ pixels strictly inside the deep DWL interior.
        dcn_seed = (c_nf > thr) & dwl_interior;

        min_area = round(500 * scaling^2);
        dcn_seed = bwareafilt(dcn_seed > 0, [min_area Inf]);
        dcn_seed = imclose(dcn_seed, strel('disk', 5*scaling));
        dcn_seed = imfill(dcn_seed, 'holes');

        % Re-clip to the deep interior — closing/fill above CAN push pixels
        % outside dwl_interior; this guarantees DCN never touches EGL/ML/IGL.
        dcn_seed = dcn_seed & dwl_interior;

        % Final size filter (post-clip — small bridge fragments removed)
        dcn_final = bwareafilt(dcn_seed > 0, [round(2000*scaling^2) Inf]);
    end
end
fprintf('DCN: %d connected components, %.2f%% of image\n', ...
    max(bwlabel(dcn_final), [], 'all'), 100*sum(dcn_final(:))/numel(dcn_final));

% Use ORIGINAL cerebellum mask (before anomaly subtraction) so anomaly
% regions get classified as anomaly rather than background.
set_bin = zeros(size(all_cerebellum_orig));
set_bin(all_cerebellum_orig == 1) = 1;

set_bin(ml_final == 1) = 5;
set_bin(pc_layer_bin == 1) = 7;
set_bin(dwl_final == 1) = 6;
set_bin(dcn_final == 1) = 8;   % DCN between DWL and IGL
set_bin(c_final == 1) = 4;
set_bin(b_final == 1) = 3;
set_bin(a_final == 1) = 2;

%% EGL/IGL POSTHOC ADJUDICATION via local NeuN intensity (V2)
% In EZH2 cKO mice the ML is shrunk, so EGL and IGL physically touch with
% no separation. The pipeline gets confused at the contact zone — in one
% folium EGL pixels get labeled IGL, elsewhere IGL pixels get labeled EGL.
% Best example: 2018_05_22_s2_2_p27.
%
% Naive per-pixel NeuN thresholding FAILS at P7. Empirically, on this
% dataset, raw NeuN by layer:
%   iEGL:  median 0.45, q95 0.87   (EGL has SIGNIFICANT NeuN at P7,
%   oEGL:  median 0.36, q95 0.51    not "near zero" as in adult)
%   IGL:   median 0.79, q95 1.00
%   ML:    median 0.30, q95 0.52
% iEGL and IGL overlap between 0.40 and 0.50, so a per-pixel global
% threshold reassigns most of EGL to IGL.
%
% Solution: two-part constraint —
%   (1) only adjudicate at the EGL/IGL contact zone (within 10 um of both
%       labels), where the failures actually live
%   (2) use very strict thresholds: only fire when raw NeuN is well outside
%       the overlap zone (0.60 for "really IGL", 0.20 for "really EGL")

c_raw = mat2gray(double(c));
c_local = imgaussfilt(c_raw, 3*scaling);   % small smoothing, bridges cell gaps
% Calibrated to the actual P7 distribution on this dataset:
%   - 95% of IGL pixels are > 0.38
%   - 25% of iEGL pixels are > 0.54
%   - we want to fire on the TAIL of EGL (apparently-IGL pixels) and the
%     TAIL of IGL (apparently-EGL pixels)
neun_high = 0.65;   % top ~10% of iEGL distribution; bulk of IGL is above
neun_low  = 0.30;   % bottom ~25% of IGL; well into iEGL territory

% No contact-zone restriction (V2 had it but missed wholesale-layer-swap
% failures where there's no nearby correctly-labeled boundary). Instead use
% TOPOLOGICAL constraint via distance from the cerebellum boundary (pia):
%   - EGL territory: pia_dist < 50 um
%   - IGL territory: pia_dist > 30 um (deeper, with overlap zone 30-50 um)
%
% Adjudicate only when the current label disagrees with BOTH the NeuN signal
% AND the topological position. Less ambiguity, fewer false positives.
pia_dist_um = bwdist(~all_cerebellum_orig) * 0.5119049;   % distance in microns
in_egl_zone = pia_dist_um < 50;     % EGL is within 50 um of pia normally
in_igl_zone = pia_dist_um > 30;     % IGL is at least 30 um deep

egl_now = (set_bin == 2) | (set_bin == 3);
igl_now = (set_bin == 4);

% EGL -> IGL: pixel currently labeled EGL, NeuN high, AND deep (not at surface)
egl_to_igl = egl_now & (c_local > neun_high) & in_igl_zone & ~in_egl_zone;
set_bin(egl_to_igl) = 4;

% IGL -> EGL: pixel currently labeled IGL, NeuN low, AND near surface
igl_to_egl = igl_now & (c_local < neun_low) & in_egl_zone & ~in_igl_zone;
set_bin(igl_to_egl & (a_nf > a_level)) = 2;
set_bin(igl_to_egl & ~(a_nf > a_level)) = 3;

fprintf('EGL/IGL adjudication V3 (topology+NeuN): EGL->IGL %d px (%.3f%%), IGL->EGL %d px (%.3f%%)\n', ...
    sum(egl_to_igl(:)), 100*sum(egl_to_igl(:))/numel(set_bin), ...
    sum(igl_to_egl(:)), 100*sum(igl_to_egl(:))/numel(set_bin));

%% RIBBON DISCONTINUITY DETECTION (V4 — DIRECTIONAL, 3 RIBBONS)
% EGL, IGL, and ML are all continuous ribbons. Per the user, the typical
% EZH2-cKO failure is a CHUNK of one or part of a lobule (~100-300 um wide)
% where the ribbon "crosses over" to the wrong label.
%
% V4 V1 used DISC closing — that expands the ribbon perpendicular AND
% along its direction equally, pulling in too much adjacent tissue.
% V4 V2 uses MULTI-ANGLE LINE STREL closing — bridges gaps ALONG the
% ribbon direction without expanding perpendicular. Per user: "we should
% expand in the predicted direction of the ribbon".
%
% Also extended:
%   - EGL reclaims from IGL, DWL, AND DCN (per user: "DWL being in the EGL"
%     and "DCN are expanding into the adjacent EGL")
%   - IGL reclaims from EGL and DWL
%   - ML now adjudicated as the third ribbon (NEW)

% pia_dist_um was computed in V3 above
line_len_um = 50;   % bridges gaps up to ~50 um along the ribbon
line_len_px = round(line_len_um / 0.5119049);

% Helper: multi-angle line-strel close — bridges gaps in any direction
% (6 line orientations every 30 deg) but doesn't expand the ribbon laterally
% the way a disc would. UNION across orientations preserves the ribbon's
% native shape while bridging breaks along its direction.
ribbon_dir_close = @(mask) ...
    imclose(mask, strel('line', line_len_px, 0))   | ...
    imclose(mask, strel('line', line_len_px, 30))  | ...
    imclose(mask, strel('line', line_len_px, 60))  | ...
    imclose(mask, strel('line', line_len_px, 90))  | ...
    imclose(mask, strel('line', line_len_px, 120)) | ...
    imclose(mask, strel('line', line_len_px, 150));

% --- EGL ribbon ---
egl_ribbon = (set_bin == 2) | (set_bin == 3);
egl_closed = ribbon_dir_close(egl_ribbon);
egl_gap = egl_closed & ~egl_ribbon & all_cerebellum_orig;
% Reclaim from IGL, DWL, OR DCN (any non-EGL deep label that landed in the
% EGL ribbon trajectory). Constraints: low NeuN AND near pia.
not_egl_target = (set_bin == 4) | (set_bin == 6) | (set_bin == 8);
egl_iEGL_gain = egl_gap & not_egl_target & (c_local < 0.35) ...
                        & (pia_dist_um < 50) & (a_nf > a_level);
egl_oEGL_gain = egl_gap & not_egl_target & (c_local < 0.35) ...
                        & (pia_dist_um < 50) & ~(a_nf > a_level);
set_bin(egl_iEGL_gain) = 2;
set_bin(egl_oEGL_gain) = 3;

% --- IGL ribbon ---
igl_ribbon = (set_bin == 4);   % refresh after EGL adjudication
igl_closed = ribbon_dir_close(igl_ribbon);
igl_gap = igl_closed & ~igl_ribbon & all_cerebellum_orig;
% Reclaim from EGL or DWL. Constraints: high NeuN AND deep tissue.
not_igl_target = (set_bin == 2) | (set_bin == 3) | (set_bin == 6);
igl_gain = igl_gap & not_igl_target & (c_local > 0.55) ...
                   & (pia_dist_um > 50);
set_bin(igl_gain) = 4;

% --- ML ribbon (NEW) ---
ml_ribbon = (set_bin == 5);
ml_closed = ribbon_dir_close(ml_ribbon);
ml_gap = ml_closed & ~ml_ribbon & all_cerebellum_orig;
% ML signal: low p27 + low-moderate NeuN + sits between EGL and IGL.
% Reclaim from IGL primarily (most common confusion), with constraint:
% NeuN moderate (not high IGL), p27 low (not EGL), pia_dist intermediate.
ml_gain = ml_gap & (set_bin == 4) & (c_local < 0.40) ...
                 & ~(a_nf > a_level) & (pia_dist_um > 30);
set_bin(ml_gain) = 5;

% --- DWL fingers (NEW) ---
% DWL is structurally a TREE/branching pattern (the cerebellar white matter
% arbor), not a smooth ribbon — but the same multi-angle directional close
% bridges along the local branch direction. The user observed: "the tips of
% the fingers are getting misassigned" — closing extends each finger along
% its trajectory and reclaims any ML tips that should be DWL.
%
% DWL signal: dim across all 3 channels (low DAPI - sparse cells, no p27,
% no NeuN). Use raw-channel combined intensity for discrimination — ML is
% brighter (parallel fibers + interneurons keep DAPI/NeuN moderate).
combined_raw = (mat2gray(double(a)) + mat2gray(double(b)) + mat2gray(double(c))) / 3;

dwl_ribbon = (set_bin == 6);
dwl_closed = ribbon_dir_close(dwl_ribbon);
dwl_gap = dwl_closed & ~dwl_ribbon & all_cerebellum_orig;
% Reclaim from ML where combined intensity is low (DWL-like) and deep
dwl_gain = dwl_gap & (set_bin == 5) & (combined_raw < 0.25) & (pia_dist_um > 50);
set_bin(dwl_gain) = 6;

fprintf('Ribbon discontinuity (V4 directional): EGL+%d (%.3f%%), IGL+%d (%.3f%%), ML+%d (%.3f%%), DWL+%d (%.3f%%)\n', ...
    sum(egl_iEGL_gain(:)) + sum(egl_oEGL_gain(:)), ...
    100*(sum(egl_iEGL_gain(:)) + sum(egl_oEGL_gain(:)))/numel(set_bin), ...
    sum(igl_gain(:)), 100*sum(igl_gain(:))/numel(set_bin), ...
    sum(ml_gain(:)), 100*sum(ml_gain(:))/numel(set_bin), ...
    sum(dwl_gain(:)), 100*sum(dwl_gain(:))/numel(set_bin));

%% DCN ANATOMICAL CONSTRAINT
% Deep Cerebellar Nuclei are by definition DEEP — they live in the white
% matter at the base of the cerebellum. A DCN label within ~50 um of the
% pial surface is anatomically impossible. These are almost certainly
% DCN-detector over-extension (the close+fill operations expanding the
% NeuN+ DWL seed past the dwl_interior boundary into adjacent EGL).
%
% The ribbon adjudication can't catch this because DCN's high NeuN signal
% trips the "EGL must have low NeuN" gate. Use topology alone.
dcn_at_surface = (set_bin == 8) & (pia_dist_um < 50);
% Reassign to iEGL or oEGL based on local p27
set_bin(dcn_at_surface & (a_nf > a_level)) = 2;
set_bin(dcn_at_surface & ~(a_nf > a_level)) = 3;
fprintf('DCN at surface reclaim: %d px (%.3f%%)\n', ...
    sum(dcn_at_surface(:)), 100*sum(dcn_at_surface(:))/numel(set_bin));

%% DWL ANATOMICAL CONSTRAINT
% DWL is the deep white matter — it should never appear at the surface.
% User reported: "2018_05_22_s2_5 shows persistent DWL on the outer surface".
% Same anatomical impossibility as DCN-at-surface: DWL within ~30 um of
% pia is detector over-extension into EGL territory.
% (DWL gets a tighter threshold than DCN because DWL legitimately sits
% closer to surface in some folium tips than DCN does.)
dwl_at_surface = (set_bin == 6) & (pia_dist_um < 30);
% Reassign to iEGL/oEGL based on local p27 — these are surface pixels
set_bin(dwl_at_surface & (a_nf > a_level)) = 2;
set_bin(dwl_at_surface & ~(a_nf > a_level)) = 3;
fprintf('DWL at surface reclaim: %d px (%.3f%%)\n', ...
    sum(dwl_at_surface(:)), 100*sum(dwl_at_surface(:))/numel(set_bin));

%% LABEL-1 RECLAIM (missing-value propagation, depth-gated)
% Pixels still labeled "1" (cerebellum, no specific layer) didn't get any
% layer assignment. Assign them via NEAREST labeled-layer neighbor — but
% gate by ANATOMICAL DEPTH so we don't propagate EGL into deep DWL or
% propagate DWL into EGL territory.
%
% User: "we are now putting an EGL across the base of the cerebellum where
% there is DWL". Original V1 (committed in 552db4f) was depth-blind —
% nearest neighbor could pick EGL even for deep pixels.
%
% V2 depth-gated rule:
%   pixel pia_dist < 30 um (surface)  -> any of {iEGL, oEGL, IGL}
%   pixel pia_dist 30-80 um (mid)     -> any of {iEGL, oEGL, ML, IGL, PCL}
%   pixel pia_dist > 80 um (deep)     -> only {IGL, ML, DWL, PCL, DCN} (NO EGL)
% PCL and DCN allowed at all depths (they're cell-body classes that can
% appear within the relevant ribbon regardless of position).
unassigned = (set_bin == 1);
n_unassigned = sum(unassigned(:));
if n_unassigned > 0
    % Allowed-label sets per depth zone
    surface_labels  = [2, 3, 4, 7];          % iEGL, oEGL, IGL, PCL
    mid_labels      = [2, 3, 4, 5, 6, 7, 8]; % all biological
    deep_labels     = [4, 5, 6, 7, 8];       % IGL, ML, DWL, PCL, DCN (NO EGL)

    % Compute zone masks
    surface_zone = pia_dist_um <  30;
    mid_zone     = (pia_dist_um >= 30) & (pia_dist_um <= 80);
    deep_zone    = pia_dist_um >  80;

    % For each zone, compute nearest allowed label
    nearest_dist  = inf(size(set_bin));
    nearest_label = zeros(size(set_bin), 'uint8');

    % Surface zone candidates
    nd_s = inf(size(set_bin));  nl_s = zeros(size(set_bin), 'uint8');
    for L = surface_labels
        layer_mask = (set_bin == L);
        if any(layer_mask(:))
            d = bwdist(layer_mask);
            closer = d < nd_s;
            nd_s(closer) = d(closer);
            nl_s(closer) = L;
        end
    end
    % Mid zone candidates
    nd_m = inf(size(set_bin));  nl_m = zeros(size(set_bin), 'uint8');
    for L = mid_labels
        layer_mask = (set_bin == L);
        if any(layer_mask(:))
            d = bwdist(layer_mask);
            closer = d < nd_m;
            nd_m(closer) = d(closer);
            nl_m(closer) = L;
        end
    end
    % Deep zone candidates (no EGL allowed)
    nd_d = inf(size(set_bin));  nl_d = zeros(size(set_bin), 'uint8');
    for L = deep_labels
        layer_mask = (set_bin == L);
        if any(layer_mask(:))
            d = bwdist(layer_mask);
            closer = d < nd_d;
            nd_d(closer) = d(closer);
            nl_d(closer) = L;
        end
    end

    % Combine: pick the right nearest_label based on zone
    nearest_label(surface_zone) = nl_s(surface_zone);
    nearest_dist(surface_zone)  = nd_s(surface_zone);
    nearest_label(mid_zone)     = nl_m(mid_zone);
    nearest_dist(mid_zone)      = nd_m(mid_zone);
    nearest_label(deep_zone)    = nl_d(deep_zone);
    nearest_dist(deep_zone)     = nd_d(deep_zone);

    % Distance limit: only assign within 50 um of a valid neighbor
    max_dist_px = round(50 / 0.5119049);    % ~98 px
    can_assign = unassigned & (nearest_dist <= max_dist_px);
    set_bin(can_assign) = nearest_label(can_assign);
    fprintf('Label-1 reclaim (depth-gated): %d unassigned, %d propagated (%.3f%%)\n', ...
        n_unassigned, sum(can_assign(:)), 100*sum(can_assign(:))/numel(set_bin));
end

%% ORPHAN FRAGMENT CLEANUP (continuity-driven)
% Continuity tests revealed ML/IGL/EGL break into many small components
% (largest CC often only 15-50% of total). Many are "orphans" — small
% isolated chunks of the wrong label inside a different layer's territory.
%
% For each ribbon layer (EGL union, IGL, ML), find small isolated components
% (< 2000 px^2 = ~520 um^2) and reassign each to its dominant boundary
% neighbor label. This addresses the ML-in-DWL and IGL-in-EGL fragments.

% (orphan_cleanup is defined as a local function at end of file)
orphan_max_area = round(500 / (0.5119049^2));    % ~1909 px (500 um^2)

n_before = nnz(set_bin == 5);
set_bin = orphan_cleanup(set_bin, [5], orphan_max_area);
ml_orphans_reassigned = n_before - nnz(set_bin == 5);

n_before = nnz(set_bin == 4);
set_bin = orphan_cleanup(set_bin, [4], orphan_max_area);
igl_orphans_reassigned = n_before - nnz(set_bin == 4);

n_before = nnz(set_bin == 2 | set_bin == 3);
set_bin = orphan_cleanup(set_bin, [2 3], orphan_max_area);
egl_orphans_reassigned = n_before - nnz(set_bin == 2 | set_bin == 3);

fprintf('Orphan cleanup: ML %d, IGL %d, EGL %d (each is total px reassigned)\n', ...
    ml_orphans_reassigned, igl_orphans_reassigned, egl_orphans_reassigned);

set_bin(pc_bin_filt == 1) = 7;
% Anomaly regions OVERRIDE all biological labels — they're untrustworthy
set_bin(anomaly_mask == 1) = 9;
%set_bin(pia_outline == 1) = 0;
figure
imagesc(set_bin)

col_map =  [1 1 1
            0.9333333 0.1725490 0.1725490
            0.0000000 0.8039216 0.4000000
            0.9333333 0.5098039 0.3843137
            1 1 1
            1 1 1
            0.3607843 0.6745098 0.9333333
            1 0 1
            0.3 0.3 0.3];   % NEW: anomaly — dark gray
show_bin = set_bin; 
b_inc = imdilate(b_final, strel('disk', 8*scaling));
show_bin(b_inc == 1) = 3;      
show_bin(ml_final == 1) = 5;
show_bin(a_final == 1) = 2;
show_bin(pc_bin_filt == 1) = 7;
show_bin(pc_layer_bin == 1) = 7;

imshow(label2rgb(show_bin, col_map, 'w'))            



% Save raw labels FIRST so we always have them even if label2rgb fails
imwrite(uint8(set_bin), [file(1:(end-4)) '_labels.tif'])
% Use explicit 10-color jet palette to handle DCN (8) and anomaly (9)
imwrite(label2rgb(set_bin, jet(10), 'k'),[file(1:(end-4)) '_segments.tif'])

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
region_names = {'all_cerebellum', 'iEGL', 'oEGL', 'IGL', 'ML', 'DWL', 'PCL', 'DCN', 'anomaly'};
% Truncate or pad to match table height (some slides may have fewer classes)
nrows = height(T);
if nrows <= numel(region_names)
    T.region_names = region_names(1:nrows)';
else
    rn_padded = region_names;
    for k = (numel(region_names)+1):nrows; rn_padded{k} = sprintf('extra_%d', k); end
    T.region_names = rn_padded(1:nrows)';
end

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
region_names = {'all_cerebellum', 'iEGL', 'oEGL', 'IGL', 'ML', 'DWL', 'PCL', 'DCN', 'anomaly'};
% Truncate or pad to match table height (some slides may have fewer classes)
nrows = height(T);
if nrows <= numel(region_names)
    T.region_names = region_names(1:nrows)';
else
    rn_padded = region_names;
    for k = (numel(region_names)+1):nrows; rn_padded{k} = sprintf('extra_%d', k); end
    T.region_names = rn_padded(1:nrows)';
end

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


function newlabels = orphan_cleanup(set_bin_in, ribbon_labels, max_area_px)
% Returns updated set_bin with orphan components of these labels reassigned.
% ribbon_labels: vector of labels that together form one "ribbon"
%                (e.g., [2 3] for iEGL+oEGL = EGL)
% max_area_px:   components <= this size are "orphans" candidate for relabel
newlabels = set_bin_in;
ribbon_mask = false(size(set_bin_in));
for L = ribbon_labels
    ribbon_mask = ribbon_mask | (set_bin_in == L);
end
cc = bwconncomp(ribbon_mask);
if cc.NumObjects == 0; return; end
stats = regionprops(cc, 'Area');
areas = [stats.Area];
orphan_idx = find(areas <= max_area_px);
if isempty(orphan_idx); return; end

for k = orphan_idx
    pixels = cc.PixelIdxList{k};
    % Build boundary mask: 1-pixel-wide ring around the orphan
    orphan_mask = false(size(set_bin_in));
    orphan_mask(pixels) = true;
    boundary = imdilate(orphan_mask, strel('disk', 2)) & ~orphan_mask;
    % Count boundary labels (exclude 0 background and same-ribbon labels)
    b_labels = newlabels(boundary);
    b_labels = b_labels(b_labels > 0);
    for L = ribbon_labels
        b_labels = b_labels(b_labels ~= L);
    end
    if isempty(b_labels); continue; end
    % Mode: most common boundary label
    new_label = mode(double(b_labels));
    newlabels(pixels) = new_label;
end
end


