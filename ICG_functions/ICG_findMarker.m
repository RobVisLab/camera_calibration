function [marker_corners, marker_similarity] = ICG_findMarker (image, template, corner_threshold, min_area, intensity_threshold)
% [marker_corners, marker_similarity] = ICG_findMarker (image, template, corner_threshold, min_area, intensity_threshold)

if nargin < 5,
    intensity_threshold = [];
end;

if size(image,3) ~= 1
    error('Image should only have one plane')
end

if size(template) ~= [32 32]
    error('Template picture dimenstion must be 32x32!')
end

[corners_col, corners_row] = findMarkerCorners(image, corner_threshold, min_area, intensity_threshold);

if isempty(corners_col),
    marker_corners = [];
    marker_similarity = 0;
    return;
end;

warning('off', 'Images:maketform:conditionNumberofAIsHigh');
matches = GetOrientationIdentity(image, template, corners_col, corners_row);
warning('on', 'Images:maketform:conditionNumberofAIsHigh');
[marker_corners, marker_similarity] = selectBestMatch(matches);
marker_corners = marker_corners';
if marker_similarity < .7
    % refine marker corners and calculate similarity
    % [marker_corners, marker_similarity] = ICG_refineMarkerCorners(image, marker_corners, template);
end
