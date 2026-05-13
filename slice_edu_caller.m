%   =======================================================================================
%   Copyright (C) 2013  Erlend Hodneland
%   Email: erlend.hodneland@biomed.uib.no 
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   =======================================================================================


% load the data
%load /Users/jpurzner/Dropbox/imaging_analysis/cellsegm/data/nucleistain_3D.mat
clear prm

% for visualization
plane = 8;

%
% Segmentation by iterative thresholding
%
% No smoothing to save time
prm.smoothim.method = 'none';

% No ridge filtering
prm.filterridges = 0;

% No illumination
prm.illum = 0;

prm.segmct.method = 'thrs';

prm.thrs.th = 0.5; 



% ridge filtering
prm.segmct.filterridges = 0;

% smoothing method
%prm.smoothim.method = 'eed';

% minima method
%prm.segmct.getminima.method = 'nucleus';

% threshold in minima
%prm.segmct.getminima.nucleus.segmct.thrs.th = 0.2;

% split cells in minima
%prm.segmct.getminima.nucleus.segmct.split = 1;

% threshold for splitting cells in minima
%prm.segmct.getminima.nucleus.segmct.splitth = 1;




% subtract nucleus channel from surface stain
imsegm = double(z);
%prm.method = 'thrs';

%prm.smoothim.eed.deltat = 1;
%prm.smoothim.eed.kappa = 10;

%prm.thrs.th = 0.5;
%imsm = cellsegm.smoothim(imsegm,'eed','prm',prm);
[cellbw1,wat,imsegmout,prmout] = cellsegm.segmct(imsegm,0.01,10,'prm',prm);

cellsegm.show(imsegm(:,:,plane),1);
title('Hoechst staining');axis off;
cellsegm.show(cellbw1(:,:,plane),2);
title('Cell segmentation, iterative thresholding, without splitting');axis off;

% Add splitting of cells. Can do this separately to have better control
% prm.splitth = 0.5;

% cells above this threshold are split (all cells here)
% n = prmout.minvolvox;
%h = [0.5 0.5 1.5];
%cellbw2 = cellsegm.splitcells(cellbw1,prm.splitth,n,h);
%cellsegm.show(cellbw2(:,:,plane),3);
%title('Cell segmentation, iterative thresholding, with splitting');axis off;
