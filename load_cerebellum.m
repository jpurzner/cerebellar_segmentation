data = bfopen('/Users/jpurzner/Dropbox/images/IF_ACK_p27/2017_09_23_ezh2_edu_4.zvi');
series20 = data{20, 1};
t = {series20{1:12, 1}};
t = cellfun(@(x) imadjust(x,stretchlim(x),[]), t, 'UniformOutput', false);
y = reshape([t{:}], [size(t{1},1), size(t{1},2), length(t)]);
z = mat2gray(y);
%J = arrafun(@x imadjust(y(:,:,5),stretchlim(y(:,:,5)),[]);