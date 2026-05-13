function egl = cerebellum_egl(dapi)  
I = rgb2gray(dapi_test);
bin = imbinarize(I, 'adaptive');
MN   =fix([6 6]);
se = strel('rectangle', MN);

bin = imdilate(bin,strel('disk',round(1),0));
BWco = imopen(bin, se);
%BWc0 =imfill(BWco,'holes');
BWco = imerode(BWco,strel('disk',round(1),0));
BWco = imopen(BWco, se);
BWco = bwareafilt(BWco,20);
%BWc  = imclose(BWc, se);
%BWco = imopen(BWc, se);
BWl = bwlabel(BWco);
CC = bwconncomp(BWco);

props = regionprops(BWl,'ConvexArea');  
idx = ( [props.ConvexArea] > 3000);
BWf = ismember(BWl,find(idx));
BWfl = bwlabel(BWf);
coloredLabels = label2rgb (BWfl, 'hsv', 'k', 'shuffle');

figure, imshowpair(coloredLabels, I, 'montage')
figure, imshow(imfuse(coloredLabels, I,'blend'))

end