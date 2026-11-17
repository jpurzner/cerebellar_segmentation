I = imread('coins.png');

[gx,gy] = gradient(double(I));
G = sqrt(gx.*gx+gy.*gy);
G = G/max(max(G));

sc = 3; % scale
nor = 8; % norientations
rad = 28; % radius

A = acc(G,rad,sc,nor);

threshold = 0.5;

centers = [];
C = zeros(size(A));

LocMax = imregionalmax(A);
[lm_rows, lm_cols] = find(LocMax == 1);
for i = 1:numel(lm_rows)
    if A(lm_rows(i),lm_cols(i)) > threshold
        centers = [centers; [lm_rows(i) lm_cols(i)]];
        C(lm_rows(i),lm_cols(i)) = 1;
    end
end

for n = 1:size(centers,1)
    for alpha = linspace(0,2*pi,100)
        x = centers(n,1)+rad*cos(alpha);
        y = centers(n,2)+rad*sin(alpha);
        C(round(x),round(y)) = 1;
    end
end

imshow([G A C])