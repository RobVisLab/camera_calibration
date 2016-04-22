function wimg = quadtobox_h(xgrid, ygrid, img, dst, M, tmplt_x_grid, tmplt_y_grid, f_type)
% QUADTOBOX_HOMO - Warp contents of a quadrilateral to a rectangle
%   W = QUADTOBOX_HOMO(IMG, DST, M, F_TYPE)
%
%   Destination corners, DST = [x1, x2, x3, x4; y1, y2, y3, y4];
%   These are the (assumed rectangular) corner points of the template
%   image and must be integers.
%
%   The projection matrix, M transforms DST into the image IMG 
%   [u' v' w] (image) = M * [x y 1] (DST - template);
%   u = u' / w; v = v' / w;
%
%   Matlab style indices, i.e. start from 1.

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: quadtobox_h.m,v 1.5 2008-09-19 12:45:39 ruether Exp $

% Check args
if nargin<8 f_type = 'bilinear'; end
if nargin<7 error('Invalid args'); end

% Dimensions of destination image - integers, assume rectangle
minv = min(dst');
maxv = max(dst');

% Get all points in destination to sample
% [tmplt_x_grid tmplt_y_grid] = meshgrid(1:maxv(1), 1:maxv(2));

xy = [reshape(tmplt_x_grid, prod(size(tmplt_x_grid)), 1)'; reshape(tmplt_y_grid, prod(size(tmplt_y_grid)), 1)'];
xy = [xy; ones(1,size(xy,2))];

% Transform into source
uv = M * xy;

% Divide for homography
uv = uv ./ repmat(uv(3,:),3,1);

% Remove homogeneous
uv = uv(1:2,:)';

% Sample
xi = reshape(uv(:,1),maxv(2),maxv(1));
yi = reshape(uv(:,2),maxv(2),maxv(1));
wimg = ICG_interp2(xgrid, ygrid, img, xi, yi, f_type);

% Check for NaN background pixels - replace them with a background of 0
idx = find(isnan(wimg));
if ~isempty(idx)
	wimg(idx) = 0;
end
