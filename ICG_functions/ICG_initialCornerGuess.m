function [corner_guess, grid_coords] = ICG_initialCornerGuess (image, marker_corners, target_type, feature_x_grid, feature_y_grid, verbose);

% ICG_initialCornerGuess Pixel accurate guess of corner positions.
%
% [corner_guess, grid_coords] = ICG_initialCornerGuess (image, marker_corners);
%
% Parameters:
%   image: image of calibration target
%   marker_corners: corner points of central marker
%
% Returns:
%   corner_guess: pixel accurate corner estimates
%   grid_coords: corner coordinates on grid
%
% Example:
%   [corner_guess, grid_coords] = ICG_initialCornerGuess (image, marker_corners);
%
% Guess of checkerboard corner positions based on correlating a pattern
% over the rectified image.
%
% Author:           Matthias R�ther
% Created:          29.12.2005
% Version:          $Revision: 1.22 $
% Matlab version:   7.0 (Win32)
%
% Used: user, user, ... (enter your name here, if you are _not_ the author and you have used 
%                       this function successfully)
%
% Copyright (c) 2005 ICG, Graz University of Technology, Austria.

if ~exist('verbose', 'var'),
    verbose = 0;
end;

if size(marker_corners, 2) ~= 4,
    error('ICG_initialCornerGuess: invalid number of marker corners.');
end;
if ~exist('target_type', 'var'),
    target_type = 'checkerboard';
end;
if ~exist('feature_x_grid', 'var') || ~exist('feature_y_grid', 'var'),
    [feature_x_grid, feature_y_grid] = meshgrid (-100:100, -100:100);
end;


[origin_offset_row, origin_offset_col] = find(feature_x_grid==0 & feature_y_grid==0);
origin_offset_xy = [origin_offset_col origin_offset_row]';
ref_feature_distance = 25;

[num_features_row, num_features_col] = size(feature_x_grid);
feature_mask_valid = im2bw(zeros(num_features_row, num_features_col));
feature_mask_failed = feature_mask_valid;
strel_3x3 = strel('square', 3);
homography_set{1}{1} = [];
image_points_x = zeros(num_features_row, num_features_col);
image_points_y = zeros(num_features_row, num_features_col);


[num_rows, num_cols] = size(image);

%% Initial Homography from Marker corners
target_corners = [-1.5 -1.5  1.5  1.5;...
                  -3.5 -0.5 -0.5 -3.5]*ref_feature_distance;
if strcmp(target_type, 'circles'),
target_corners = [-1.75 -1.75  1.75  1.75;...
                  -3.75 -0.25 -0.25 -3.75]*ref_feature_distance;   
end;

marker_corners(3,:) = 1;
target_corners(3,:) = 1;
initial_homography = ICG_fitAffinity2D(marker_corners, target_corners); % Skrabal data only works with affinity
feature_mask_failed((-3:-1)+origin_offset_xy(2), (-1:1)+origin_offset_xy(1)) = 1;

%% Initial Point Hypotheses
initial_points_grd = [-2 -1 0 1 2  -2  2 -2  2 -2  2 -2 -1  0  1  2       -3:3         -3:3 ones(1,5)*-3 ones(1,5)*3;...
                       0  0 0 0 0  -1 -1 -2 -2 -3 -3 -4 -4 -4 -4 -4  ones(1,7) ones(1,7)*-5         -4:0         -4:0];
point_indices = [initial_points_grd(1,:)+origin_offset_xy(1);initial_points_grd(2,:)+origin_offset_xy(2)];

feature_mask_rect = [1 1 size(feature_mask_valid, 2), size(feature_mask_valid, 1)];
pts_outside_grd = find(~insideRectangle(point_indices, feature_mask_rect));
initial_points_grd(:, pts_outside_grd) = [];
point_indices     (:, pts_outside_grd) = [];

initial_points_ref = initial_points_grd*ref_feature_distance;
initial_points_ref(3,:) = 1;
initial_points_img = ICG_normalizePoints(inv(initial_homography) * initial_points_ref);

image_rect = [1 1 size(image, 2), size(image, 1)];
pts_outside_img = find(~insideRectangle(initial_points_img, image_rect));
initial_points_grd(:, pts_outside_img) = [];
point_indices(:, pts_outside_img) = [];
initial_points_ref(:, pts_outside_img) = [];
initial_points_img(:, pts_outside_img) = [];

for i=1:size(point_indices, 2),
    feature_mask_valid(point_indices(2,i), point_indices(1,i)) = 1;
    image_points_x(point_indices(2,i), point_indices(1,i)) = initial_points_img(1,i);
    image_points_y(point_indices(2,i), point_indices(1,i)) = initial_points_img(2,i);
    homography_set{point_indices(2,i)}{point_indices(1,i)} = initial_homography;
end;

[image_points_x, image_points_y, feature_mask_valid, feature_mask_failed] = refinePoints(point_indices(:,17:end), ...
    image_points_x, image_points_y, homography_set, ...
    feature_x_grid, feature_y_grid, ...
    feature_mask_valid, feature_mask_failed, ...
    image, ref_feature_distance, target_type);

%% Main Loop
num_points_visited_old = 0;
num_points_visited_new = 1;
while num_points_visited_new - num_points_visited_old >= 1,
    num_points_visited_old = num_points_visited_new;
    
    %% Find candidates neighboring valid corners
    candidate_mask = xor(imdilate(feature_mask_valid, strel_3x3), feature_mask_valid) & ~feature_mask_failed;
    [point_indices_row, point_indices_col] = find(candidate_mask);
    point_indices = [point_indices_col point_indices_row]';
    
    if isempty(point_indices),
        break;
    end;  
    
    [homography_set, image_points_x, image_points_y, feature_mask_valid, feature_mask_failed] = findHomographies(point_indices,...
        image_points_x, image_points_y, homography_set, ...
        feature_x_grid, feature_y_grid, ...
        feature_mask_valid, feature_mask_failed, ...
        ref_feature_distance);
    
    [image_points_x, image_points_y, feature_mask_valid, feature_mask_failed] = refinePoints(point_indices, ...
        image_points_x, image_points_y, homography_set, ...
        feature_x_grid, feature_y_grid, ...
        feature_mask_valid, feature_mask_failed, ...
        image, ref_feature_distance, target_type);
    
    num_points_visited_new = sum(sum(feature_mask_valid | feature_mask_failed));
    
    if verbose >= 1,
        fprintf(1, '.');
    end;
end;

%% Eliminate central points
initial_points_grd = [-2 -1 0 1 2  -2  2 -2  2 -2  2 -2 -1  0  1  2  ;...
                       0  0 0 0 0  -1 -1 -2 -2 -3 -3 -4 -4 -4 -4 -4  ];
point_indices = [initial_points_grd(1,:)+origin_offset_xy(1);initial_points_grd(2,:)+origin_offset_xy(2)];

feature_mask_rect = [1 1 size(feature_mask_valid, 2), size(feature_mask_valid, 1)];
pts_outside_grd = find(~insideRectangle(point_indices, feature_mask_rect));
point_indices     (:, pts_outside_grd) = [];
for i=1:size(point_indices, 2),
    feature_mask_valid(point_indices(2,i), point_indices(1,i)) = 0;
end;


%% Compile Output
point_indices = find(feature_mask_valid);
corner_guess = [image_points_x(point_indices), image_points_y(point_indices)]';
grid_coords = [feature_x_grid(point_indices), feature_y_grid(point_indices)]';
grid_coords(3,:) = 0;


%% Find Homographies
function [homography_set, image_points_x, image_points_y, feature_mask_valid, feature_mask_failed] = findHomographies(point_indices,...
        image_points_x, image_points_y, homography_set, ...
        feature_x_grid, feature_y_grid, ...
        feature_mask_valid, feature_mask_failed, ...
        ref_feature_distance);

feature_mask_rect = [1 1 size(feature_mask_valid, 2), size(feature_mask_valid, 1)];
[nhood_5x5_xg, nhood_5x5_yg] = meshgrid(-2:2, -2:2);
new_feature_mask_valid = feature_mask_valid;
for i=1:size(point_indices, 2),
    current_point = point_indices(:,i);
    
    nhood_x = current_point(1)+nhood_5x5_xg(:)';
    nhood_y = current_point(2)+nhood_5x5_yg(:)';
    pts_outside = ~insideRectangle([nhood_x;nhood_y], feature_mask_rect);
    nhood_x(pts_outside) = [];
    nhood_y(pts_outside) = [];
    
    nhood_inside_ind = sub2ind(size(feature_mask_valid), nhood_y, nhood_x);
    valid_nhood_ind = nhood_inside_ind(feature_mask_valid(nhood_inside_ind));
    if length(valid_nhood_ind) < 4,
        failed_nhood_ind = nhood_inside_ind(feature_mask_failed(nhood_inside_ind));
        if length(valid_nhood_ind) + length(failed_nhood_ind) + sum(pts_outside) >= 24,
            feature_mask_failed(current_point(2), current_point(1)) = 1;
            new_feature_mask_valid(current_point(2), current_point(1)) = 0;
        end;
        continue;
    end;

%     valid_nhood_x = nhood_x(valid_nhood_ind);
%     valid_nhood_y = nhood_y(valid_nhood_ind);
% 
%     valid_nhood_ind = sub2ind(size(feature_mask_valid), valid_nhood_y, valid_nhood_x);
    
    corners_ref = [feature_x_grid(valid_nhood_ind);...
                   feature_y_grid(valid_nhood_ind)] * ref_feature_distance;
    corners_ref(3,:) = 1;  
    corners_img = [image_points_x(valid_nhood_ind);...
                   image_points_y(valid_nhood_ind)];
    corners_img(3,:) = 1;
    homography = ICG_fitAffinity2D(corners_img, corners_ref);
    %% FIXX Validity check of homography
    %% FIXX check det(homography)
%     if(det(homography)>10e-10)
        current_point_ref = [feature_x_grid(current_point(2), current_point(1))*ref_feature_distance;...
                             feature_y_grid(current_point(2), current_point(1))*ref_feature_distance;...
                             1];

        current_point_img = ICG_normalizePoints(inv(homography)* current_point_ref);

        new_feature_mask_valid(current_point(2), current_point(1)) = 1;
        homography_set{current_point(2)}{current_point(1)} = homography;
        image_points_x(current_point(2), current_point(1)) = current_point_img(1);
        image_points_y(current_point(2), current_point(1)) = current_point_img(2); 
%     else
%         continue;
%     end
end;
feature_mask_valid = new_feature_mask_valid;
stop=1;

%% Refine Point
function [image_points_x, image_points_y, feature_mask_valid, feature_mask_failed] = refinePoints(point_indices, ...
        image_points_x, image_points_y, homography_set, ...
        feature_x_grid, feature_y_grid, ...
        feature_mask_valid, feature_mask_failed, ...
        image, ref_feature_distance, target_type);

min_similarity = 0.7;
min_intensity_difference = 0.05;
feature_mask_grid = [feature_x_grid(1,1) feature_y_grid(1,1) size(feature_mask_valid, 2) size(feature_mask_valid, 1)];
feature_mask_img = [1 1 size(image, 2) size(image, 1)];
  
    
persistent refinePoints_template;
persistent refinePoints_target_type;
persistent refinePoints_ref_feature_distance;
if isempty(refinePoints_template) || ~strcmp(refinePoints_target_type, target_type) || refinePoints_ref_feature_distance ~= ref_feature_distance,
    refinePoints_target_type = target_type;
    refinePoints_ref_feature_distance = ref_feature_distance;
    refinePoints_template = createCorrTemplate (target_type, ref_feature_distance);
end;

% persistent image_xg;
% persistent image_yg;
% [xgr, xgc] = size(image_xg);
% [ygr, ygc] = size(image_yg);
% [r, c] = size(image);
% if isempty(image_xg) || isempty(image_yg) || xgr ~= r || xgc ~= c || ygr ~= r || ygc ~= c,
%     [image_xg image_yg] = meshgrid(size(image));
% end;


min_neighbor_area_image = 10*10;

global USE_IPP;
if isempty(USE_IPP),
    USE_IPP = 0;
end;

if USE_IPP,
    image_uint8_t = im2uint8(image)';
%     uint8(image*255)';
end;

for i=1:size(point_indices, 2),
    current_point = point_indices(:,i);
    if ~feature_mask_valid(current_point(2), current_point(1)),
        continue;
    end;
    homography = homography_set{current_point(2)}{current_point(1)};

    current_point_grid = [feature_x_grid(current_point(2), current_point(1));...
                          feature_y_grid(current_point(2), current_point(1))];
    neighbor_points_grid = [[-1 -1  1 1]+current_point_grid(1); ...
                            [ 1 -1 -1 1]+current_point_grid(2)];
                        
    pts_outside_grd = ~insideRectangle(neighbor_points_grid, feature_mask_grid);
    if any(pts_outside_grd),
        feature_mask_failed(current_point(2), current_point(1)) = 0;
        feature_mask_valid(current_point(2), current_point(1)) = 1;
        continue;
    end;
                            
    neighbor_points_ref(1:2,:) = neighbor_points_grid(1:2,:) * ref_feature_distance;
    neighbor_points_ref(3,:) = 1;
    
    neighbor_points_img = ICG_normalizePoints(inv(homography) * neighbor_points_ref);
    pts_outside_img = ~insideRectangle(neighbor_points_img, feature_mask_img);
    if any(pts_outside_img),
        feature_mask_failed(current_point(2), current_point(1)) = 1;
        feature_mask_valid(current_point(2), current_point(1)) = 0;
        continue;
    end;
    
    neighbor_area_image = polyarea(neighbor_points_img(1,:), neighbor_points_img(2,:));
    
    if neighbor_area_image < min_neighbor_area_image,
        feature_mask_failed(current_point(2), current_point(1)) = 1;
        feature_mask_valid(current_point(2), current_point(1)) = 0;
        continue;
    end;
    
    if USE_IPP,
        src_roi = [1 1 size(image,2) size(image,1)];
        dst_roi = [min(neighbor_points_ref(1,:)) min(neighbor_points_ref(2,:)) max(neighbor_points_ref(1,:)) max(neighbor_points_ref(2,:))];
        dst_roi(3:4) = dst_roi(3:4) - dst_roi(1:2) + 1;
        
        translation = eye(3);
        translation(1:2,3) = [1-dst_roi(1) 1-dst_roi(2)];
        dst_roi(1:2) = [1 1];
        
        [warped_image status] = ICG_projectiveWarpIPP(image_uint8_t, translation*homography, src_roi, dst_roi, 1);
        
        if status ~= 0,
            tform = maketform('projective', homography');
            warped_image = imtransform(image, tform, 'XData', [min(neighbor_points_ref(1,:)) max(neighbor_points_ref(1,:))], ...
                                                     'YData', [max(neighbor_points_ref(2,:)) min(neighbor_points_ref(2,:))]);
        end;
            
        warped_image = flipdim(warped_image, 1);
        warped_image = im2double(warped_image);
        
%         tform = maketform('projective', homography');
%         warped_image1 = imtransform(image, tform, 'XData', [min(neighbor_points_ref(1,:)) max(neighbor_points_ref(1,:))], ...
%                                                  'YData', [max(neighbor_points_ref(2,:)) min(neighbor_points_ref(2,:))]);
% 
%         d = warped_image1 - warped_image;
%         if any(abs(d(:))>1.1/255),
%             stop = 1;
%         end;
    else,
        tform = maketform('projective', homography');
        warped_image = imtransform(image, tform, 'XData', [min(neighbor_points_ref(1,:)) max(neighbor_points_ref(1,:))], ...
                                                 'YData', [max(neighbor_points_ref(2,:)) min(neighbor_points_ref(2,:))]);
    end;
                                             
    
    [offset_cr, similarity, contrast] = calcOffset (refinePoints_template, warped_image);
    if contrast < min_intensity_difference,
        feature_mask_failed(current_point(2), current_point(1)) = 1;
        feature_mask_valid(current_point(2), current_point(1)) = 0;
        continue;
    end;
    if similarity < min_similarity,
        feature_mask_failed(current_point(2), current_point(1)) = 1;
        feature_mask_valid(current_point(2), current_point(1)) = 0;
        continue;
    end;
    current_point_ref = current_point_grid * ref_feature_distance;
    refined_point_ref = current_point_ref;
    refined_point_ref(1) = current_point_ref(1)+offset_cr(1);
    refined_point_ref(2) = current_point_ref(2)-offset_cr(2);
    refined_point_ref(3) = 1;
    refined_point_img = ICG_normalizePoints(inv(homography)*refined_point_ref);
    
    image_points_x(current_point(2), current_point(1)) = refined_point_img(1);
    image_points_y(current_point(2), current_point(1)) = refined_point_img(2);
    feature_mask_valid(current_point(2), current_point(1)) = 1;
end;


%--------------------------------------------------------------------------
function [offset_cr, similarity, contrast] = calcOffset (template, search_window);

[t_r, t_c] = size(template);
[s_r, s_c] = size(search_window);

corr_image = abs(ICG_normXCorr2(template, search_window, 'valid'));
[c_r, c_c] = size(corr_image);
diff_r = s_r - c_r;
diff_c = s_c - c_c;
corr_image = padarray(corr_image, [floor(diff_r/2), floor(diff_c/2)], 0.0, 'both');

% start_index_r = 1+floor(t_r/2);
% end_index_r = c_r-floor(t_r/2);
% start_index_c = 1+floor(t_c/2);
% end_index_c = c_c-floor(t_c/2);
% corr_image = corr_image (start_index_r:end_index_r, start_index_c:end_index_c);
[c_r, c_c] = size(corr_image);

corr_image = abs(corr_image);
[num_rows, num_cols] = size (corr_image);
maximum = max(corr_image(:));
if maximum == 0,
    offset_cr = [];
    similarity = 0;
    contrast = 0;
    return;
end;
[max_pos_r, max_pos_c] = find(corr_image == maximum);
max_pos_r = max_pos_r(1);
max_pos_c = max_pos_c(1);
offset_cr = [max_pos_c, max_pos_r] - [ceil(c_c/2), ceil(c_r/2)];
similarity = maximum;
region_r = floor(max_pos_r-t_r/2):floor(max_pos_r-t_r/2)+t_r;
region_c = floor(max_pos_c-t_c/2):floor(max_pos_c-t_c/2)+t_c;
if region_r(1) == 0,
    region_r = region_r+1;
end;
if region_c(1) == 0,
    region_c = region_c+1;
end;
found_region = search_window(region_r, region_c);
contrast = max(found_region(:)) - min(found_region(:));
offset_cr = offset_cr';

%--------------------------------------------------------------------------
function template = createCorrTemplate (template_type, ref_feature_distance);
% template = createCorrTemplate (template_type, ref_feature_distance);

template_size = round(2*ref_feature_distance * 0.75);
if mod(template_size, 2) ~= 0,
    template_size = template_size-1;
end;
half_template_size = template_size/2;

if strcmp(template_type, 'checkerboard'),
    template = zeros(template_size,template_size);
    template(1:half_template_size, 1:half_template_size) = 1;
    template(half_template_size+1:end, half_template_size+1:end) = 1;
elseif strcmp(template_type, 'circles'), 
    template = ones(template_size,template_size);
    disk = fspecial('disk', ref_feature_distance/4);
    disk = disk/max(disk(:));
    disk = abs(disk-1);
    template(1:size(disk,1), 1:size(disk,2)) = disk;
    template(end+1-size(disk,1):end, 1:size(disk,2)) = disk;
    template(1:size(disk,1), end+1-size(disk,2):end) = disk;
    template(end+1-size(disk,1):end, end+1-size(disk,2):end) = disk;
end;

%--------------------------------------------------------------------------
function is_inside = insideRectangle(points, rectangle_xywh);

is_inside = points(1,:)>=rectangle_xywh(1) & points(1,:)<=rectangle_xywh(1)-1+rectangle_xywh(3) & ...
            points(2,:)>=rectangle_xywh(2) & points(2,:)<=rectangle_xywh(2)-1+rectangle_xywh(4);
        