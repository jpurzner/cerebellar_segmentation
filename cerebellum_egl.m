function egl = cerebellum_egl(dapi)  



Is = imread(dapi);
%Is = rgb2gray(I);
I = imgaussfilt(Is,4);
h = fspecial('log', [20 20], 0.3);
I = imfilter(I,h,'replicate');
%test = edge(I,'prewitt', 0.3);
%figure, imshow(test)
bin = imbinarize(I, 'adaptive', 'Sensitivity',0.5);
MN   =fix([6 6]);
se = strel('rectangle', MN);



BWd = imdilate(bin,strel('disk',round(1),0));
BWco = imopen(BWd, se);
%BWco =imfill(BWco,'holes');
BWco = imerode(BWco,strel('disk',round(1),0));
BWco = imopen(BWco, se);
BWco = bwareafilt(BWco,20);
%BWc  = imclose(BWc, se);
%BWco = imopen(BWc, se);
BWl = bwlabel(BWco);
CC = bwconncomp(BWco);

props = regionprops(BWl,'ConvexArea');  
idx = ( [props.ConvexArea] > 2000);
BWf = ismember(BWl,find(idx));
BWfl = bwlabel(BWf);

% minimize convex hull area by eliminating 

bw = activecontour(I,BWf, 50, 'edge');

coloredLabels = label2rgb (BWfl, 'hsv', 'k', 'shuffle');

figure, imshowpair(imfuse(coloredLabels, Is,'blend'), I, 'montage')
figure, imshowpair(bin, bw, 'montage')

egl = BWfl;


    function area = hull_area(bw) 
        CH = bwconvhull(BW)
        
    end % end hull_area 

end % end cerebellum_egl