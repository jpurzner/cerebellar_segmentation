function layer = cerebellum_layers(file)  

a = imread(file,1);
b = imread(file, 2);
c = imread(file, 3);

ha = homomorphic(a, 2, 0.25,2,0,5);
hb = homomorphic(b, 2, 0.25,2,0,5);
hc = homomorphic(c, 2, 0.25,2,0,5);

ha_g =  imgaussfilt(ha,5);
hb_g =  imgaussfilt(hb,5);
hc_g =  imgaussfilt(hc,5);

col_g(:,:,1) = ha_g;
col_g(:,:,2) = hb_g;
col_g(:,:,3) = hc_g;

[L,N] = superpixels(col_g,20000);


a_n = a;
a_n(a_n< 8000) =0;

b_n = b;
b_n(b_n< 10000) = NaN;


c_n = c;
c_n(c_n< 10000) = NaN;


a_norm = homomorphic(a_n, 2, 0.25,2,0,5);
%a_norm = (double(a_n)./double(b_n)) + double(a_n) ;
%a_norm(a_norm == Inf) = NaN;
%a_norm(a_norm == -Inf) = NaN;

c_norm = (double(c_n)./double(b_n)) + (2*double(c_n)) ;
c_norm(c_norm == Inf) = NaN;
c_norm(c_norm == -Inf) = NaN;

dapi_high_mask = b; 
dapi_high_mask(b >= 42500) = 0;
dapi_high_mask(b < 42500) = 1;



%hb_n = homomorphic(b_n, 2, 0.25,2,0,5);
%ha = homomorphic(a, 2, 0.25,2,0,5);
b_norm = (double(b_n)./double(a_n))   ;
%b_norm = (double(hb_n)./double(ha))   ;
b_norm(b_norm == Inf) = NaN;
b_norm(b_norm == -Inf) = NaN;

%layer = b_norm;
%return

b_norm(isnan(b_norm))=0;
b_norm = b_norm.^2;
b_norm = ordfilt2(b_norm,16,ones(4,4));
b_norm =  imgaussfilt(b_norm,4);
b_norm = illumination_cor(b_norm);
b_norm = mat2gray(b_norm);

%layer = b_norm;
%return

b_bin = b_norm > 0.05;

b_bin = imerode(b_bin,strel('disk',round(2),0));
b_bin = imdilate(b_bin,strel('disk',round(2),0));
b_bin = bwareaopen(b_bin, 2000);
b_bin = imclearborder(b_bin, 4);

%s = regionprops(b_bin,'area');




c_norm = c_norm .* (1 - double(b_bin));

c_norm(isnan(c_norm))=0;
c_norm = c_norm.^2;
c_norm = ordfilt2(c_norm,16,ones(4,4));
c_norm = imgaussfilt(c_norm,5);
c_norm = mat2gray(c_norm);
c_bin = c_norm > 0.2;
%diff_im = anisodiff2D(c_norm_c, 15,1/7,30,2);
c_bin = imdilate(c_bin,strel('disk',round(2),0));
c_bin = imerode(c_bin,strel('disk',round(2),0));
c_bin = imdilate(c_bin,strel('disk',round(4),0));
c_bin = imerode(c_bin,strel('disk',round(4),0));
%c_bin = imdilate(c_bin,strel('disk',round(6),0));
%c_bin = imerode(c_bin,strel('disk',round(6),0));
c_bin = bwmorph(c_bin,'thick', 8);
c_bin = imdilate(c_bin,strel('disk',round(4),0));
c_bin = imerode(c_bin,strel('disk',round(4),0));
c_bin = bwmorph(c_bin,'thick', 8);
c_bin = imerode(c_bin,strel('disk',round(8),0));
c_bin = bwmorph(c_bin,'thick', 10);
c_bin = imerode(c_bin,strel('disk',round(10),0));
c_bin = bwmorph(c_bin,'thick', 12);
c_bin = imerode(c_bin,strel('disk',round(12),0));
c_bin = imdilate(c_bin,strel('disk',round(10),0));
c_bin = imerode(c_bin,strel('disk',round(10),0));

%c_bin = imdilate(c_bin,strel('disk',round(8),0));
%c_bin = imerode(c_bin,strel('disk',round(8),0));

%c_bin = imdilate(c_bin,strel('disk',round(8),0));
%c_bin = imerode(c_bin,strel('disk',round(8),0));
%c_bin =imfill(c_bin,'holes');
%c_bin = imdilate(c_bin,strel('disk',round(10),0));
%c_bin = imerode(c_bin,strel('disk',round(10),0));



c_bin = bwareaopen(c_bin, 50000);
c_bin = imclearborder(c_bin, 4);
c_bin = imfill(c_bin, 'holes');

%s = regionprops(c_bin,'area');




bb = mat2gray(b);
all_cerebellum = bb > 0.1;
all_cerebellum = imclearborder(all_cerebellum, 8);
all_cerebellum = imfill(all_cerebellum, 'holes');


% seperate the deep white matter from the molecular layer 
ha_n = homomorphic(a_n, 2, 0.25,2,0,5);
% this seem trivial but you can't have too many zeros in the denominator 
% or the ratio image looks weird
hc = homomorphic(c, 2, 0.25,2,0,5);
dwl_norm = log2(ha_n./hc);
dwl_norm = anisodiff(dwl_norm, 15, 50, 0.1, 1) ;



dwl_norm = mat2gray(dwl_norm);
dwl_norm = dwl_norm .* (1-c_bin);
dwl_norm = dwl_norm.* (1-b_bin);
% use a very high cutoff so can later seperate the 


dwl_bin = dwl_norm >= 0.85;
dwl_bin = imerode(dwl_bin,strel('disk',round(8),0));
dwl_bin = imdilate(dwl_bin,strel('disk',round(2),0));
dwl_bin = imerode(dwl_bin,strel('disk',round(2),0));
dwl_bin = imdilate(dwl_bin,strel('disk',round(4),0));
dwl_bin = imerode(dwl_bin,strel('disk',round(4),0));
%c_bin = imdilate(c_bin,strel('disk',round(6),0));
%c_bin = imerode(c_bin,strel('disk',round(6),0));
dwl_bin = bwmorph(dwl_bin,'thick', 8);
dwl_bin = imdilate(dwl_bin,strel('disk',round(4),0));
dwl_bin = imerode(dwl_bin,strel('disk',round(4),0));
dwl_bin = bwmorph(dwl_bin,'thick', 8);
dwl_bin = imerode(dwl_bin,strel('disk',round(8),0));
%dwl_bin = bwmorph(dwl_bin,'thick', 10);
%dwl_bin = imdilate(dwl_bin,strel('disk',round(5),0));
dwl_bin = imfill(dwl_bin, 'holes');



[L, num] = bwlabel(dwl_bin, 8);
count_pixels_per_obj = sum(bsxfun(@eq,L(:),1:num));
[~,ind] = max(count_pixels_per_obj);
dwl_bin = (L==ind);
dwl_bin = imdilate(dwl_bin,strel('disk',round(4),0));
dwl_bin = imfill(dwl_bin, 'holes');
%layer = dwl_bin;
%return

 

% get the molecular layer
hb = homomorphic(b, 2, 0.25,2);
ml = hb;



%ml = (b/1000).^2;
ml = anisodiff(ml, 15, 75, 0.1, 2);
%ml =  imgaussfilt(ml,5);
ml =  mat2gray(ml);
ml = 1 - ml;
ml = ml .* all_cerebellum;
ml = ml .* (1-c_bin);
ml = ml .* (1-b_bin);
ml = ml .* (1-dwl_bin);

ml  = ml > 0.7;

ml = imdilate(ml,strel('disk',round(2),0));
ml = imerode(ml,strel('disk',round(2),0));
ml = imdilate(ml,strel('disk',round(4),0));
ml = imerode(ml,strel('disk',round(4),0));
%c_bin = imdilate(c_bin,strel('disk',round(6),0));
%c_bin = imerode(c_bin,strel('disk',round(6),0));
ml = bwmorph(ml,'thick', 8);
ml = imdilate(ml,strel('disk',round(4),0));
ml = imerode(ml,strel('disk',round(4),0));
ml = bwmorph(ml,'thick', 8);
ml = imerode(ml,strel('disk',round(8),0));
%dwl_bin = bwmorph(dwl_bin,'thick', 10);
%dwl_bin = imdilate(dwl_bin,strel('disk',round(5),0));
ml = imfill(ml, 'holes');



%ml =  mat2gray(ml);
%ml = 1 - ml;
%ml = ml .* all_cerebellum;
%ml = ml .* (1-c_bin);
%ml = ml .* (1-b_bin);
%ml =  imgaussfilt(ml,5);




ml_lab = bwlabel(ml);

%
%ml_lab(c_bin)  = (max(max(ml_lab)) + 1);
%ml_lab((a_bin | b_bin) & (~ ml_lab) )  = (max(max(ml_lab)) + 1);
%ml_lab(b_bin  & (~ ml_lab) )  = (max(max(ml_lab)) + 1);

%layer = ml_lab;
%return


% crudely segregate iEGL
%a_norm = a_norm .* (1 - double(c_bin));
a_norm = anisodiff(a_norm, 15, 75, 0.1, 2);

a_norm(isnan(a_norm))=0;


%a_norm = a_norm.^2;
%a_norm = ordfilt2(a_norm,16,ones(4,4));
%a_norm = imgaussfilt(a_norm,5);
a_norm = mat2gray(a_norm);

a_bin = a_norm  > 0.75;
a_bin = a_bin .* (1 - double(c_bin));
a_bin = a_bin .* (1 - double(ml));
a_bin = a_bin .* (1 - double(dwl_bin));
a_bin = a_bin == 1;

ml_lab = bwlabel(ml);
ml_lab(c_bin)  = (max(max(ml_lab)) + 1);
%ml_lab((a_bin | b_bin) & (~ ml_lab) )  = (max(max(ml_lab)) + 1);
ml_lab(b_bin)  = (max(max(ml_lab)) + 1);
ml_lab(dwl_bin)  = (max(max(ml_lab)) + 1);
ml_lab(a_bin)  = (max(max(ml_lab)) + 1);

all_lab = double(ml); 
all_lab(dwl_bin) = 2;
all_lab(c_bin)  = 3;
all_lab(b_bin)  = 4;
all_lab(a_bin) = 5;

col_n(:,:,1) = mat2gray(a_norm, [0.6 0.9]);
col_n(:,:,2) = mat2gray(b_norm, [0 0.4]);
col_n(:,:,3) = mat2gray(c_norm, [0.1 0.4]);

layer = all_lab;
return

a_bin = a_norm > 0.1;
a_bin = bwareaopen(a_bin, 5000);
a_bin = imclearborder(a_bin, 4);
a_bin = imerode(a_bin,strel('disk',round(1),0));
a_bin = imdilate(a_bin,strel('disk',round(1),0));
%c_bin = imdilate(c_bin,strel('disk',round(6),0));
%c_bin = imerode(c_bin,strel('disk',round(6),0));
a_bin = bwmorph(a_bin,'thick', 2);
a_bin = imerode(a_bin,strel('disk',round(2),0));
a_bin = bwareaopen(a_bin, 5000);


% check to see if any of the regions in a_bin are misassigned
a_lab = bwlabel(a_bin);
a_lab(c_bin)  = (max(max(a_lab)) + 1);





ml = bwareaopen(ml, 5000);
ml = imclearborder(ml, 4);
ml = imerode(ml,strel('disk',round(4),0));
ml = imdilate(ml,strel('disk',round(4),0));
%c_bin = imdilate(c_bin,strel('disk',round(6),0));
%c_bin = imerode(c_bin,strel('disk',round(6),0));
ml = bwmorph(ml,'thick', 4);
ml = imerode(ml,strel('disk',round(4),0));
ml = bwareaopen(ml, 5000);









c_label = all_cerebellum;
c_label = c_label + b_bin;
c_label(a_bin) = 3;
c_label(c_bin) = 4;





% no gap between the oEGL and the iEGL

% must be at least 1 pixel iEGL between molecular layer





layer = c_label;
return




end 
