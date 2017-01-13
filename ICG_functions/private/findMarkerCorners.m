function [corners_col, corners_row] = findMarkerCorners(image, corner_threshold, min_area, intensity_threshold)

    image = im2double(image);
    image = image./max(image(:));
    
    if nargin < 4 || isempty(intensity_threshold),
        pattern_bwimage = adaptivethreshold(image, min_area/10);
    else
        pattern_bwimage = im2bw(image, intensity_threshold);
    end
    
    [connectivity_image,area] = ConnectedComponents(pattern_bwimage, min_area);
    B = bwboundaries(connectivity_image,8,'noholes');
    [B,corners] = lineSegments(B,area,corner_threshold);

    warning('off', 'MATLAB:singularMatrix');
    [corners_col corners_row] = calcIntersection(B,corners);
    warning('on', 'MATLAB:singularMatrix');
    [num_rows, num_cols] = size(image);

    bad_indices = [];
    for i=1:size(corners_col,1),
        current_region_col = corners_col(i,:);
        current_region_row = corners_row(i,:);

        % Are we outside of the image with one of the *corners*?
        if sum(current_region_col<0 | current_region_row<0 | current_region_col>num_cols | current_region_row>num_rows)>0,
            bad_indices = [bad_indices, i];
            continue;
        end
        
        % Are we suspiciously close to one of the image borders with any *edge*?
        edge_too_close_threshold = 7;
        if sum(current_region_row - edge_too_close_threshold < 0) > 1 || ...
           sum(current_region_row + edge_too_close_threshold > num_rows) > 1 || ...
           sum(current_region_col - edge_too_close_threshold < 0) > 1 || ...
           sum(current_region_col + edge_too_close_threshold > num_cols) > 1
           bad_indices = [bad_indices, i];
           %disp(sprintf('####### dropping marker %d, too close to the border', i));
        end
    end
    corners_col(bad_indices,:) = [];
    corners_row(bad_indices,:) = [];
