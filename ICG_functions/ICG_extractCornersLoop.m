function [is_good, corner_points, grid_coords, verbose_info] = ICG_extractCornersLoop (image, template, parameters, verbose_info, index);
% [is_good, corner_points, grid_coords, verbose_info] = ICG_extractCornersLoop (image, template, parameters, verbose_info, index);

global USE_IPP;
backup_USE_IPP = USE_IPP;
USE_IPP = parameters.USE_IPP;

min_marker_similarity = 0.4;
current_image = im2double(image);

if parameters.verbose >= 1,
    fprintf(1, 'Finding Marker...');
end;
% if(numel(parameters.approx_marker_width_pixels)>1) % we want to search an area
marker_similarity = 0;
marker_corners = [];
for i=1:numel(parameters.approx_marker_width_pixels)
    [curr_marker_corners, curr_marker_similarity] = ICG_findMarker (current_image, template, ...
    parameters.approx_marker_width_pixels(i)*2/3, parameters.approx_marker_width_pixels(i)*parameters.approx_marker_width_pixels(i), parameters.image_threshold);
    if(curr_marker_similarity>marker_similarity)
        marker_similarity=curr_marker_similarity;
        marker_corners = curr_marker_corners;
    end
end
% else
%     [marker_corners, marker_similarity] = ICG_findMarker (current_image, template, ...
%         parameters.approx_marker_width_pixels*2/3, parameters.approx_marker_width_pixels*parameters.approx_marker_width_pixels, parameters.image_threshold);
%     % marker_corners = marker_corners'; 
% end

if parameters.verbose >= 1,
    fprintf(1, 'DONE. Similarity: %f\n', marker_similarity);
end;

if marker_similarity < min_marker_similarity,
    is_good = 0;
    corner_points = [];
    grid_coords = [];
    USE_IPP = backup_USE_IPP;
    return;
end;

if parameters.verbose >= 1,
    fprintf(1, 'Initial corner guess');
end;
[feature_x_grid, feature_y_grid] = meshgrid(parameters.grid_coordinates_h, parameters.grid_coordinates_v);
[corner_guess, grid_coord] = ICG_initialCornerGuess (current_image, marker_corners, parameters.target_type, feature_x_grid, feature_y_grid, parameters.verbose);
if isempty(corner_guess),
    num_corners = 0;
else,
    num_corners = size(corner_guess, 2);
end;  
if parameters.verbose >= 1,
    fprintf(1, 'DONE. Found %d corners.\n', num_corners);
end;


%% DEBUG
% if ~isempty(corner_guess),
%     figure, imshow(current_image);
%     hold on, plot(corner_guess(1,:), corner_guess(2,:), 'r.');
% end;
% is_good = 0;
% corner_points = [];
% grid_coords = [];
% verbose_info = [];
% USE_IPP = backup_USE_IPP;
% return;
%% NODEBUG

if parameters.verbose >= 1,
    fprintf(1, 'Corner refinement');
end;
[corners, grid_coordinates, checker_verbose_info] = ICG_checkerboardCorners (current_image, corner_guess, grid_coord, ...
    parameters.corner_extraction_method, parameters.target_type, parameters.verbose, parameters.USE_IPP, parameters.feature_subsampling_rate);
if isempty(corner_guess),
    num_corners = 0;
else,
    num_corners = size(corners, 2);
end;  
if parameters.verbose >= 1,
    fprintf(1, 'DONE. %d corners remain.\n', num_corners);
end;

verbose_info.figure_handles{index} = checker_verbose_info.figure_handles;
if isfield(parameters, 'verbose') && parameters.verbose >= 2,
    drawnow;
    pause(0.2);
end;
if ~isempty(corners),
    corners(3,:) = 1;
end;

if size(corners, 2) < 6,
    if isfield(parameters, 'verbose') && parameters.verbose >= 2,
        handle = figure; 
        imshow(current_image); title('Bad Image');
        verbose_info.figure_handles{index} = handle; 
    end;
    is_good = 0;
    corner_points = [];
    grid_coords = [];
    USE_IPP = backup_USE_IPP;
    return;
end;
is_good = 1;
corner_points = corners;
grid_coords = grid_coordinates * parameters.grid_width_mm;
