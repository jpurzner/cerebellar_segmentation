
%


outputFileName = 'img_stack2.tif';
delete img_stack2.tif;
for K=1:length(wat(1, 1, :))
   coloredLabels = label2rgb (wat(:,:,K), 'hsv', 'k', 'shuffle');
   imwrite(coloredLabels, outputFileName, 'WriteMode', 'append',  'Compression','none');
   %imwrite(wat(:, :, K), outputFileName, 'WriteMode', 'append',  'Compression','none');
end

outputFileName = 'img_stack3.tif';
delete img_stack3.tif
for K=1:length(wat(1, 1, :))
   extractmask=bwmorph(wat(:,:,K) > 0,'remove');
   RGB=y(:,:,K);
   RGB(:,:,2)=extractmask;
   RGB(:,:,3)=1;
   imwrite(RGB, outputFileName, 'WriteMode', 'append',  'Compression','none');
   %imwrite(wat(:, :, K), outputFileName, 'WriteMode', 'append',  'Compression','none');
end