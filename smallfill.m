function  new  = smallfill( original, sz )
%smallfill removes small holes 
% from https://blogs.mathworks.com/steve/2008/08/05/filling-small-holes/
filled = imfill(original, 'holes');
holes = filled & ~original;
bigholes = bwareaopen(holes, sz);
smallholes = holes & ~bigholes;
new = original | smallholes;

end

