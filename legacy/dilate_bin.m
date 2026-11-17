function bin = dilate_bin(im_bin, sz, rem)
    im_bin = imclose(im_bin, strel('disk',sz));
    im_bin = imfill(im_bin, 'holes');
    im_bin = imdilate(im_bin, strel('disk',sz));
    im_bin = bwareafilt(im_bin > 0,[rem Inf]);
    bin = im_bin; 
    return 
end 