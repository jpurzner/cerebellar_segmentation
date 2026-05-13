function new_im = background_subtract(im) 
        [pixelCounts , grayLevels] = histcounts(im,500);
        [max_pixelCounts, max_i] =  max(pixelCounts(2:end));
        % get high point of histogram (noise) and add the typical width 
        low_in = min(grayLevels(max_i + 1) + 0.04, 0.9)  ;
        % set saturation where 95% of non zero pixels
        % clip the low pixel counts 
        low_i = min(find(grayLevels >= low_in));
        high_i = min(find(cumsum(pixelCounts(low_i:end))/sum(pixelCounts(low_i:end)) > 0.9));
        high_i = low_i + high_i ;
        high_in = min([grayLevels(high_i) 1]);
        new_im = imadjust(im, [low_in high_in], [0 1]);
        return
end