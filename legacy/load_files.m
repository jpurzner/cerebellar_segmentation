a = imread('/Users/jpurzner/Dropbox/images/2017_05_18_p27_p7_ezh2_cko/2017_05_17_C_p27_dapi_10x_pano_red_cor_crop.tif');
b = imread('/Users/jpurzner/Dropbox/images/2017_05_18_p27_p7_ezh2_cko/2017_05_17_C_p27_dapi_10x_pano_red_cor_crop.tif', 2);
c = imread('/Users/jpurzner/Dropbox/images/2017_05_18_p27_p7_ezh2_cko/2017_05_17_C_p27_dapi_10x_pano_red_cor_crop.tif', 3);
ha = homomorphic(a, 2, 0.25,2,0,5);
hb = homomorphic(b, 2, 0.25,2,0,5);
hc = homomorphic(c, 2, 0.25,2,0,5);

ha_g =  imgaussfilt(ha,5);
hb_g =  imgaussfilt(hb,5);
hc_g =  imgaussfilt(hc,5);

oegl = hb_g./ha_g;
oegl = mat2gray(oegl, [0.6 2]); 
iegl = mat2gray(ha_g, [0.4 0.7]);
igl = mat2gray(hc_g, [0.4 0.7]);

clear col_g
col_g(:,:,1) = iegl;
col_g(:,:,2) = oegl;
col_g(:,:,3) = igl;

clear col
col(:,:,1) = mat2gray(ha_g, [0.4 0.7]);
col(:,:,2) = iegl;
col(:,:,3) = igl;
