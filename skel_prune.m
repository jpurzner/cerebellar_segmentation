function skelD = skel_prune(ima) 

% Teja Muppirala
% https://www.mathworks.com/matlabcentral/answers/88284-remove-the-spurious-edge-of-skeleton
%

skel= bwmorph(ima,'skel',Inf);
B = bwmorph(skel, 'branchpoints');
E = bwmorph(skel, 'endpoints');
[y,x] = find(E);
B_loc = find(B);
Dmask = false(size(skel));
for k = 1:numel(x)
    D = bwdistgeodesic(skel,x(k),y(k));
    distanceToBranchPt = min(D(B_loc));
    Dmask(D < distanceToBranchPt) =true;
end
skelD = skel - Dmask;
return

end 