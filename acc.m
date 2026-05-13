% Computes the accumulator image for a
% Wavelet-Based Circular Hough Transform, given a particular radius
%
% REFERENCE
% Marcelo C., Davi G., and Kris G.
% Wavelet-Based Circular Hough Transform
% and its Application in Embryo Development Analysis.
% VISAPP, Barcelona, 2013. BibTeX:
% @proceedings{Cicconet2013,
%   editor    = {Sebastiano Battiato and
%                Jos{\'e} Braz},
%   title     = {VISAPP 2013 - Proceedings of the International Conference
%                on Computer Vision Theory and Applications, Volume 1, Barcelona,
%                Spain, 21-24 February, 2013},
%   booktitle = {VISAPP (1)},
%   publisher = {SciTePress},
%   year      = {2013},
%   isbn      = {978-989-8565-47-1}
% }
%
% USAGE
% A = acc(G,rad,sc,nor)
%
% RETURNS
% A: the accumulator image (in the range [0,1])
%
% PARAMETERS
% G: gradient image (see example below)
% rad: radius of circles we're looking for
% sc: scale of the wavelet filter to be used;
%    this parameter should be approximately
%    the the size of the edges of the circles in the image
% nor: number of orientations of the wavelet to be used;
%    typical values are 8 and 16
%
% VERSION
% 1.0, Mar 21 2013
%
% AUTHOR
% Marcelo Cicconet
% marceloc.net

function A = acc(G,rad,sc,nor)

[nr,nc] = size(G);

rrs = zeros(2,nor);
crs = zeros(2,nor);

rr1 = rad+1;
rr2 = nr-rad;
cr1 = rad+1;
cr2 = nc-rad;

S = zeros(rr2-rr1+1,cr2-cr1+1);

for or = 1:nor
    ang = (or-1)/nor*2*pi;
    rrs(1,or) = rr1+round(rad*cos(ang));
    rrs(2,or) = rr2+round(rad*cos(ang));
    crs(1,or) = cr1+round(rad*sin(ang));
    crs(2,or) = cr2+round(rad*sin(ang));
end

Js = zeros(size(S,1),size(S,2),nor);
for or = 1:nor
    J = G(rrs(1,or):rrs(2,or),...
        crs(1,or):crs(2,or));
    
    ang = (or-1)/nor*360;
    [mr,~] = morlet(sc,ang,1);

    J = conv2(J,mr,'same');
    Js(:,:,or) = J.*(J > 0);
    
    S = S+J;
end

A = zeros(nr,nc);
S = (S/max(max(S))).^2;
A(rr1:rr2,cr1:cr2) = S;

end