function normalized_points = ICG_normalizePoints( points );

% ICG_normalizePoints Normalize points such that the homogeneous part is 1.
%
% normalized_points = ICG_normalizePoints( points );
%
% Parameters:
%   points: set of n points of dimension m in homogeneous representation, format
%           mxn matrix.
%
% Returns:
%   normalized_points: set of n points of dimension m in homogeneous
%                      representation, format mxn matrix. The last row of the
%                      matrix is 1.
%
% Example:
%   points2d = [1 2;...
%               2 4;...
%               1 2];
%   norm_points2d = ICG_normalizePoints( points2d );
%
% Normalize a set of points in homogeneous representation, such that the
% homogeneous part of each point is 1.
%
% See also:
%
% Author:           Matthias Rüther
% Created:          12.07.2002
% Version:          $Revision: 1.3 $
% Matlab version:   6.1
%
% Used:
%
% Copyright (c) 2002 ICG, Graz University of Technology, Austria.

% $RCSfile: ICG_normalizePoints.m,v $ $Revision: 1.3 $ $Date: 2006-01-19 10:27:30 $
%
% Revision Information:
%
%

[num_coords, num_points] = size(points);

% If homogeneous coordinate is complex, make it real.
conj_hom_coord = conj(points(num_coords,:));
normalized_points = points .* repmat(conj_hom_coord, [num_coords,1]);

% Divide by homogeneous coordinate.
hom_coord = normalized_points(num_coords,:);
normalized_points = normalized_points ./ repmat(hom_coord, [num_coords,1]);
