s4_C_p27 = cerebellum_threshold_segment2('/Users/jpurzner/Dropbox/images/2017_05_18_p27_p7_ezh2_cko/2017_05_17_C_p27_dapi_10x_pano_red_cor_crop.tif');
s4_A_p27 = cerebellum_threshold_segment2('/Users/jpurzner/Dropbox/images/2017_05_18_p27_p7_ezh2_cko/a_p27_cor-prune.tif');
% has a tear 
s4_B_p27 = cerebellum_threshold_segment2('/Users/jpurzner/Dropbox/images/2017_05_18_p27_p7_ezh2_cko/s4_B_p27_cor_merge_crop2.tif');
s4_E_p27 = cerebellum_threshold_segment2('/Users/jpurzner/Dropbox/images/2017_05_18_p27_p7_ezh2_cko/2017_05_17_E_p27_dapi_10x_pano_cor_merge_crop.tif');
s4_F_p27 = cerebellum_threshold_segment2('/Users/jpurzner/Dropbox/images/2017_05_18_p27_p7_ezh2_cko/s4_F_p27_merge_crop.tif');
s4_G_p27 = cerebellum_threshold_segment2('/Users/jpurzner/Dropbox/images/2017_05_18_p27_p7_ezh2_cko/2017_05_17_G_p27_dapi_10x_pano_red_cor_merge_crop.tif');
close all
s5_A_p27 = cerebellum_threshold_segment2('/Users/jpurzner/Dropbox/images/2017_05_18_p27_p7_ezh2_cko/s5_A_p27_cor_merge_crop.tif');
s5_B_p27 = cerebellum_threshold_segment2('/Users/jpurzner/Dropbox/images/2017_05_18_p27_p7_ezh2_cko/s5_B_p27_cor_merge_crop.tif');
s5_C_p27 = cerebellum_threshold_segment2('/Users/jpurzner/Dropbox/images/2017_05_18_p27_p7_ezh2_cko/s5_C_p27_cor_merge_crop.tif');
s5_E_p27 = cerebellum_threshold_segment2('/Users/jpurzner/Dropbox/images/2017_05_18_p27_p7_ezh2_cko/s5_E_p27_cor_merge_crop.tif');
s5_F_p27 = cerebellum_threshold_segment2('/Users/jpurzner/Dropbox/images/2017_05_18_p27_p7_ezh2_cko/s5_F_p27_cor_merge_crop.tif');
s5_G_p27 = cerebellum_threshold_segment2('/Users/jpurzner/Dropbox/images/2017_05_18_p27_p7_ezh2_cko/s5_G_p27_cor_merge_crop.tif');
close all
%s6_B tissue disrupted
s6_C_p27 = cerebellum_threshold_segment2('/Users/jpurzner/Dropbox/images/2017_05_18_p27_p7_ezh2_cko/s6_CC_p27_merge_crop.tif');
s6_E_p27 = cerebellum_threshold_segment2('/Users/jpurzner/Dropbox/images/2017_05_18_p27_p7_ezh2_cko/s6_E_p27_cor_merge_crop.tif');
s6_G_p27 = cerebellum_threshold_segment2('/Users/jpurzner/Dropbox/images/2017_05_18_p27_p7_ezh2_cko/s6_G_p27_merge_crop.tif'); 



s4_B_p27.slide = repmat({'s4_B_p27'}, size(s4_B_p27,1),1);
s4_E_p27.slide = repmat({'s4_E_p27'}, size(s4_E_p27,1),1);
s4_F_p27.slide = repmat({'s4_F_p27'}, size(s4_F_p27,1),1);
s4_G_p27.slide = repmat({'s4_G_p27'}, size(s4_G_p27,1),1);

s5_A_p27.slide = repmat({'s5_A_p27'}, size(s5_A_p27,1),1);
s5_B_p27.slide = repmat({'s5_B_p27'}, size(s5_B_p27,1),1);
s5_C_p27.slide = repmat({'s5_C_p27'}, size(s5_C_p27,1),1);
s5_E_p27.slide = repmat({'s5_E_p27'}, size(s5_E_p27,1),1);
s5_F_p27.slide = repmat({'s5_F_p27'}, size(s5_F_p27,1),1);
s5_G_p27.slide = repmat({'s5_G_p27'}, size(s5_G_p27,1),1);

s6_C_p27.slide = repmat({'s6_C_p27'}, size(s6_C_p27,1),1);
s6_E_p27.slide = repmat({'s6_E_p27'}, size(s6_E_p27,1),1);
s6_G_p27.slide = repmat({'s6_G_p27'}, size(s6_G_p27,1),1);

cerebellum_area = [s4_B_p27; s4_E_p27; s4_F_p27; s4_G_p27; ...
      s5_A_p27; s5_B_p27; s5_C_p27; s5_E_p27; s5_F_p27; s5_G_p27; ...
      s6_C_p27; s6_E_p27; s6_G_p27 ];
  
writetable(cerebellum_area,'~/Dropbox/imaging_analysis/cerebellar_segmentation/cerebellum_area_ezh2ko.txt')  
  