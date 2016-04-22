function wimg = warp_h(xgrid, ygrid, img, p, dst, tmplt_x_grid, tmplt_y_grid)
% WARP_H - Projective warp the image
%   WIMG = WARP_H(xgrid, ygrid, IMG, P, DST)
%   Warp image IMG to WIMG. DST are the destination points, i.e. the corners
%   of the template image. P are the affine warp parameters that project
%   DST into IMG.
%
%   P = [p1, p4, p7        (i.e. transpose of what is in LK20-1)
%        p2, p5, p8   
%        p3, p6, 1]; 

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: warp_h.m,v 1.3 2006-11-07 09:42:18 ruether Exp $

if nargin<7 error('Not enough input arguments'); end

M = p;
M(1,1) = M(1,1) + 1;
M(2,2) = M(2,2) + 1;

% Use bilinear filtering to warp back to template
wimg = quadtobox_h(xgrid, ygrid, img, dst, M, tmplt_x_grid, tmplt_y_grid, 'bilinear');
