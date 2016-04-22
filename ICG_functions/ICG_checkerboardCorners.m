function [corner_coords, grid_coords, verbose_info] = ICG_checkerboardCorners (image, initial_corner_positions, grid_positions, method, template_type, verbose, USE_IPP, subsampling_rate);

% ICG_checkerboardCorners Find sub-pixel corners on checkerboard pattern.
%
% [corner_coords, grid_coords, verbose_info] = ICG_checkerboardCorners (image, initial_corner_positions, grid_positions, method, template_type, verbose, USE_IPP);
%
% Parameters:
%   image: image of checkerboard pattern
%   initial_corner_positions: estimated positions of checkerboard corners
%   grid_positions: corner coordinates on checkboard target
%   method: corner extraction method, 'linecross', 'projective' or 'saddle'
%   template_type: calibration target type ('checkerboard' or 'circles')
%   verbose: flag specifying whether to display debug info
%   USE_IPP: whether to use the IPP for a faster refinement
%
% Returns:
%   corner_coords: refined corner coordinates in the image
%   grid_coords: grid coordinates of refined ocrners
%   verbose_info: structure containing verbose information
%
% Example:
%   [corner_coords, grid_coords, verbose_info] = ICG_checkerboardCorners(image, initial_corner_positions, grid_positions, method, verbose);
%
% Refines pixel-accurate estimates of corners on a rectangular grid.
% Sub-pixel accurate edgels are searched on the four edges forming the
% corner, two lines are fitted and intersected to find the corner.  
%
% See also:
%   ICG_initialCornerGuess
%
% Author:           Matthias R�ther
% Created:          29.12.2005
% Version:          $Revision: 1.34 $
% Matlab version:   7.0
%
% Used: user, user, ... (enter your name here, if you are _not_ the author and you have used 
%                       this function successfully)
%
% Copyright (c) 2005 ICG, Graz University of Technology, Austria.

if exist('subsampling_rate', 'var')==0 || isempty(subsampling_rate),
    subsampling_rate = 1;
end;

if exist('verbose', 'var')==0 || isempty(verbose),
    verbose = 0;
end;

if exist('method', 'var')==0
    method = 'linecross';
end

if ~exist('USE_IPP','var') || isempty(USE_IPP),
    USE_IPP = 0;
end;

[num_rows, num_cols] = size(image);
num_corners = size(initial_corner_positions, 2);

neighbor_indices = zeros(1,4);

cc = ones(2,num_corners)*NaN;
gc = ones(3,num_corners)*NaN;
corner_coords = ones(2,num_corners)*NaN;
grid_coords = ones(3,num_corners)*NaN;



valid_mask_x = mod(grid_positions(1,:),subsampling_rate)==0;
valid_mask_y = mod(grid_positions(2,:),subsampling_rate)==0;
valid_mask = valid_mask_x & valid_mask_y;


%% Saddle point method
if strcmp(method, 'saddle'),
    himage = actxserver('HalconX.HImageX');
    imwrite(image, 'temp.bmp');
    himage.invoke('ReadImage', 'temp.bmp');
    [saddle_r,saddle_c] = himage.invoke('SaddlePointsSubPix', 'gauss', 1.2, 4);
    saddle_points = [saddle_c+1; saddle_r+1];
    real_image = himage.invoke('ConvertImageType', 'real');
    corner_response = real_image.invoke('CornerResponse', 3, 0.04);
    corner_response = corner_response.invoke('ScaleImageMax');
    corner_response.invoke('WriteImage', 'tiff', 0, 'corner_response.tif');
    himage.delete;
    cr = imread('corner_response.tif');
    !del temp.bmp
    !del corner_response.tif
 
    for corner_index = 1:num_corners,
         current_corner = initial_corner_positions(:,corner_index);
         current_gridpoint = grid_positions(:,corner_index);
         distances = sqrt((saddle_points(1,:)-current_corner(1)).^2 + (saddle_points(2,:)-current_corner(2)).^2);
         candidate_indices = find(distances < 5);
         if isempty(candidate_indices),
             continue;
         end;
         candidates = saddle_points(:, candidate_indices);
         corner_responses = candidates(1,:) * 0;
         for i=1:size(candidates, 2),
            corner_responses(i) = cr(round(candidates(2,i)), round(candidates(1,i)));
         end;
         [max_response mri] = max(corner_responses);
         mri = candidate_indices(mri(1));
         cc(:,corner_index) = saddle_points(:, mri);
         gc(:,corner_index) = current_gridpoint;
    end;
    
%% Line crossing method
elseif strcmp(method, 'linecross'),
    if verbose >= 1,       
       num_lines = 10;
       linelength = max(ceil(num_corners/num_lines), 100);
    end;
    for corner_index = 1:num_corners,
        if verbose>=1,
            if mod(corner_index, linelength) == 0,
                fprintf(1, '%d/%d\n', corner_index, num_corners);
            end;
            if mod(corner_index, 10) == 1,
                fprintf(1, '.');
            end;
        end;


        current_corner = initial_corner_positions(:,corner_index);
        current_gridpoint = grid_positions(:,corner_index);

        differences(1,:) = grid_positions(1,:) - current_gridpoint(1);
        differences(2,:) = grid_positions(2,:) - current_gridpoint(2);

        neighbor_indices = [];
        neighbor_indices = [neighbor_indices, find(differences(1,:)== 0 & differences(2,:)== 1)]; % upper index
        neighbor_indices = [neighbor_indices, find(differences(1,:)== 0 & differences(2,:)==-1)]; % lower index
        neighbor_indices = [neighbor_indices, find(differences(1,:)==-1 & differences(2,:)== 0)]; % left index
        neighbor_indices = [neighbor_indices, find(differences(1,:)== 1 & differences(2,:)== 0)]; % right index
        if length(neighbor_indices) ~= 4,
            continue;
        end;

        refined_corner = ICG_refineCheckerFrameLinecross (image, current_corner, initial_corner_positions(:,neighbor_indices));
        if ~isempty(refined_corner) && norm(refined_corner - current_corner) < 10,
            cc(:,corner_index) = refined_corner;
            gc(:,corner_index) = current_gridpoint;
        else,
            debug = 1;
        end;
    end;
  
%% Projective registration method
elseif strcmp(method, 'projective'),
    if verbose >= 1,       
       num_lines = 10;
       linelength = max(ceil(num_corners/num_lines), 100);
    end;

    for corner_index = 1:num_corners,
        if ~valid_mask(corner_index),
            continue;
        end;

        if verbose>=1,
            if mod(corner_index, linelength) == 0,
                fprintf(1, '%d/%d\n', corner_index, num_corners);
            end;
            if mod(corner_index, 10) == 1,
                fprintf(1, '.');
            end;
        end;

        current_corner = initial_corner_positions(:,corner_index);
        current_gridpoint = grid_positions(:,corner_index);

        differences(1,:) = grid_positions(1,:) - current_gridpoint(1);
        differences(2,:) = grid_positions(2,:) - current_gridpoint(2);
        
        neighbor_indices = [];
        neighbor_indices = [neighbor_indices, find(differences(1,:)== 0 & differences(2,:)== 1)]; % upper index
        neighbor_indices = [neighbor_indices, find(differences(1,:)== 0 & differences(2,:)==-1)]; % lower index
        neighbor_indices = [neighbor_indices, find(differences(1,:)==-1 & differences(2,:)== 0)]; % left index
        neighbor_indices = [neighbor_indices, find(differences(1,:)== 1 & differences(2,:)== 0)]; % right index
        if length(neighbor_indices) ~= 4,
            continue;
        end;

         refined_corner = ICG_refineCheckerFrameProjective (image, current_corner, initial_corner_positions(:,neighbor_indices), template_type, verbose, USE_IPP);
         if ~isempty(refined_corner) && norm(refined_corner - current_corner) < 30,
            cc(:,corner_index) = refined_corner;
            gc(:,corner_index) = current_gridpoint;
         else
            debug = 1;
        end;
    end;  
%% Projective registration method
elseif strcmp(method, 'saddle_matlab'),
    if verbose >= 1,       
       num_lines = 10;
       linelength = max(ceil(num_corners/num_lines), 100);
    end;

    for corner_index = 1:num_corners,
        if ~valid_mask(corner_index),
            continue;
        end;

        if verbose>=1,
            if mod(corner_index, linelength) == 0,
                fprintf(1, '%d/%d\n', corner_index, num_corners);
            end;
            if mod(corner_index, 10) == 1,
                fprintf(1, '.');
            end;
        end;

        current_corner = initial_corner_positions(:,corner_index);
        current_gridpoint = grid_positions(:,corner_index);

        differences(1,:) = grid_positions(1,:) - current_gridpoint(1);
        differences(2,:) = grid_positions(2,:) - current_gridpoint(2);
        
        neighbor_indices = [];
        neighbor_indices = [neighbor_indices, find(differences(1,:)== 0 & differences(2,:)== 1)]; % upper index
        neighbor_indices = [neighbor_indices, find(differences(1,:)== 0 & differences(2,:)==-1)]; % lower index
        neighbor_indices = [neighbor_indices, find(differences(1,:)==-1 & differences(2,:)== 0)]; % left index
        neighbor_indices = [neighbor_indices, find(differences(1,:)== 1 & differences(2,:)== 0)]; % right index
        if length(neighbor_indices) ~= 4,
            continue;
        end;

         refined_corner = ICG_refineCheckerFrameSaddle (image, current_corner, initial_corner_positions(:,neighbor_indices), template_type, verbose, USE_IPP);
         if ~isempty(refined_corner) && norm(refined_corner - current_corner) < 30,
            cc(:,corner_index) = refined_corner;
            gc(:,corner_index) = current_gridpoint;
         else
            debug = 1;
        end;
    end;  
%% Iterative circle correlation
elseif strcmp(method, 'circles_corr'),
    if verbose >= 1,       
       num_lines = 10;
       linelength = max(ceil(num_corners/num_lines), 100);
    end;

    for corner_index = 1:num_corners,
        if ~valid_mask(corner_index),
            continue;
        end;
        if verbose>=1,
            if mod(corner_index, linelength) == 0,
                fprintf(1, '%d/%d\n', corner_index, num_corners);
            end;
            if mod(corner_index, 10) == 1,
                fprintf(1, '.');
            end;
        end;

        current_corner = initial_corner_positions(:,corner_index);
        current_gridpoint = grid_positions(:,corner_index);

        differences(1,:) = grid_positions(1,:) - current_gridpoint(1);
        differences(2,:) = grid_positions(2,:) - current_gridpoint(2);
        
        neighbor_indices = [];
        neighbor_indices = [neighbor_indices, find(differences(1,:)== 0 & differences(2,:)== 1)]; % upper index
        neighbor_indices = [neighbor_indices, find(differences(1,:)== 0 & differences(2,:)==-1)]; % lower index
        neighbor_indices = [neighbor_indices, find(differences(1,:)==-1 & differences(2,:)== 0)]; % left index
        neighbor_indices = [neighbor_indices, find(differences(1,:)== 1 & differences(2,:)== 0)]; % right index
        if length(neighbor_indices) ~= 4,
            continue;
        end;

        refined_corner = ICG_refineCheckerFrameCirclesCorr (image, current_corner, initial_corner_positions(:,neighbor_indices), template_type, verbose, USE_IPP);
        if ~isempty(refined_corner) && norm(refined_corner - current_corner) < 30,
            cc(:,corner_index) = refined_corner;
            gc(:,corner_index) = current_gridpoint;
        else,
            debug = 1;
        end;
    end;  
%% LK registration method
elseif strcmp(method, 'LukasKanade'),  
    error('Experimental code, proved to be inaccurate. Do not use.');
    
    % Assume the target is no greater than 64x64 checkers
    [template checker_size_pix, corner_points_target, grid_points_target, grid_points_indexed_target] = ...
        ICG_createCalibrationTarget(25.4/600*64, [25.4/600*64*64 25.4/600*64*64]);
    
    [corr_points_1, corr_points_2, gp] = ICG_pointCorrespondences (...
        corner_points_target, grid_points_indexed_target, initial_corner_positions, grid_positions);
    
    img = uint8(image*255);
    img = padarray(img, (size(template) - size(img)), 'post');
    [refined_points, status, track_error] = ICG_OCV_cvCalcOpticalFlowPyrLK(uint8(~template*255), img, corr_points_1, corr_points_2, [15 15], 50, 2);

    displacements = sqrt(sum((refined_points - corr_points_2).^2, 1));
    bad = status ~= 1 | displacements > 10;
    refined_points(:, bad) = NaN;

    cc = refined_points;
    gc = grid_positions;
    gc(:, bad) = NaN;
%% Circle Fitting using Halcon lib
elseif strcmp(method, 'CircleFitting'),
    if verbose >= 1,       
       num_lines = 10;
       linelength = max(ceil(num_corners/num_lines), 100);
    end;

    
    % Estimate circle size from min distance of corners
    sqr_distances = (sum(diff(initial_corner_positions(1:2,:), 1, 2).^2, 1));
    min_distance = sqrt(min(sqr_distances)) * 0.5;
    min_circle_circumference = min_distance*0.2 *2*pi;

    % Determine circle centers
    xld = ICG_HALCON_edgesSubPix (im2uint8(image)', 'canny', 2.0, int32(10), int32(30));
    contours = ICG_HALCON_selectContoursXld (xld, {'contour_length'}, int32(min_circle_circumference), int32(99999), -0.5, 0.5);
    ICG_HALCON_clearObj(xld);
    clear xld;

    camparam = {1e6,0,1e-5,1e-5,100,100,200,200};
    [pose1, pose2] = ICG_HALCON_getCirclePose(contours, camparam, 1.0, {'center_normal'});
    ICG_HALCON_clearObj(contours);
    clear contours;
    pose1 = cat(1, pose1{:});
    pose2 = cat(1, pose2{:});
    pose1 = reshape(pose1, [6, length(pose1)/6]);
    pose2 = reshape(pose2, [6, length(pose2)/6]);
    centers1 = pose1(1:3,:);
    centers2 = pose2(1:3,:);
    normals1 = pose1(4:6,:);
    normals2 = pose2(4:6,:);
    [r1 c1] = ICG_HALCON_project3DPoint(num2cell(centers1(1,:)), num2cell(centers1(2,:)), num2cell(centers1(3,:)), camparam);
    [r2 c2] = ICG_HALCON_project3DPoint(num2cell(centers2(1,:)), num2cell(centers2(2,:)), num2cell(centers2(3,:)), camparam);
    imagepoints1 = [cat(2, c1{:}); cat(2, r1{:})];
    imagepoints1 = imagepoints1+1; 
    imagepoints1(3,:) = 1;
    imagepoints2 = [cat(2, c2{:}); cat(2, r2{:})];
    imagepoints2 = imagepoints2+1; 
    imagepoints2(3,:) = 1;
    
%     figure, quiver3(centers1(1,:), centers1(2,:), centers1(3,:), normals1(1,:), normals1(2,:), normals1(3,:));
%     figure, quiver3(centers2(1,:), centers2(2,:), centers2(3,:), normals2(1,:), normals2(2,:), normals2(3,:));
    
    imagepoints = [imagepoints1, imagepoints2];
    centers = [centers1, centers2];
    normals = [normals1, normals2];
    [plane, p, inliers] = ICG_ransacFitVector(normals, 0.05);
    inliermask = zeros(1, length(p));
    inliermask(inliers) = 1;
    mask = logical(inliermask);
    
    num_template_cols_eighth = 12;
    num_template_rows_eighth = 12;
    num_template_rows = num_template_rows_eighth*8;
    num_template_cols = num_template_cols_eighth*8;
    template = ICG_createRegistrationTemplate('circles', num_template_rows, num_template_cols, 0);
    for corner_index = 1:num_corners,
        if verbose>=1,
            if mod(corner_index, linelength) == 0,
                fprintf(1, '%d/%d\n', corner_index, num_corners);
            end;
            if mod(corner_index, 10) == 1,
                fprintf(1, '.');
            end;
        end;

        current_corner = initial_corner_positions(:,corner_index);
        current_gridpoint = grid_positions(:,corner_index);

        differences(1,:) = grid_positions(1,:) - current_gridpoint(1);
        differences(2,:) = grid_positions(2,:) - current_gridpoint(2);
        
        neighbor_indices = [find(differences(1,:)== 0 & differences(2,:)== 1),... % upper index
                            find(differences(1,:)== 0 & differences(2,:)==-1),... % lower index
                            find(differences(1,:)==-1 & differences(2,:)== 0),... % left index
                            find(differences(1,:)== 1 & differences(2,:)== 0)]; % right index
        if length(neighbor_indices) ~= 4,
            continue;
        end;
        
        refined_corner = ICG_refineCheckerFrameCircles (current_corner, initial_corner_positions(:,neighbor_indices), imagepoints(:, mask), template, verbose);
        if ~isempty(refined_corner) && norm(refined_corner - current_corner) < 10,
            cc(:,corner_index) = refined_corner;
            gc(:,corner_index) = current_gridpoint;
        else,
            debug = 1;
        end;
    end;  
end;

bad = find(isnan(cc(1,:)));
cc(:, bad) = [];
gc(:, bad) = [];
num_corners = size(cc, 2);
differences = [];

% if strcmp(method, 'linecross') || strcmp(method, 'saddle')
%     %% Impose local planarity
%     for corner_index = 1:num_corners,
%         current_corner = cc(:,corner_index);
%         current_gridpoint = gc(:,corner_index);
%          
%         differences(1,:) = gc(1,:) - current_gridpoint(1);
%         differences(2,:) = gc(2,:) - current_gridpoint(2);
%     
%         neighbor_indices = [];
%         neighbor_indices = [neighbor_indices, find(differences(1,:)== 0 & differences(2,:)== 1)]; % upper index
%         neighbor_indices = [neighbor_indices, find(differences(1,:)== 0 & differences(2,:)==-1)]; % lower index
%         neighbor_indices = [neighbor_indices, find(differences(1,:)==-1 & differences(2,:)== 0)]; % left index
%         neighbor_indices = [neighbor_indices, find(differences(1,:)== 1 & differences(2,:)== 0)]; % right index
%         if length(neighbor_indices) ~= 4,
%             continue;
%         end;
%      
%         refined_corner = ICG_imposeLocalPlanarity (image, current_corner, cc(:,neighbor_indices));
%         if ~isempty(refined_corner) && norm(refined_corner - current_corner) < 10
%             corner_coords(:,corner_index) = refined_corner;
%             grid_coords(:,corner_index) = current_gridpoint;
%         else
%             debug = 1;ICG_extractCorners
%         end
%     end
%     bad = find(isnan(corner_coords(1,:)));
%     corner_coords(:, bad) = [];
%     grid_coords(:, bad) = [];
%     
% else
    corner_coords = cc;
    grid_coords = gc;
% end

debug_image_handle = [];
verbose_info.figure_handles = [];
if verbose>=2,
    debug_image_handle = figure;
    warning('off', 'Images:truesize:imageTooBigForScreen');
    warning('off', 'Images:initSize:adjustingMag');    
    imshow(image);
    warning('on', 'Images:truesize:imageTooBigForScreen');
    warning('on', 'Images:initSize:adjustingMag');    
    figure(debug_image_handle), hold on, plot(initial_corner_positions(1,:), initial_corner_positions(2,:), 'r.');
    figure(debug_image_handle), hold on, plot(corner_coords(1,:), corner_coords(2,:), 'g.');
    if strcmp(method, 'CircleFitting'),
        figure(debug_image_handle), hold on, plot(imagepoints(1, mask), imagepoints(2, mask), 'y.');
    end;
    num_grid_points = size(grid_positions, 2);
    if verbose >=3,
        for i=1:num_grid_points,
            figure(debug_image_handle), hold on, text(initial_corner_positions(1,i)+5, initial_corner_positions(2,i)+5, sprintf('(%d|%d)', grid_positions(1,i), grid_positions(2,i)), 'FontSize',6, 'Color', 'r');
        end;
    end;
    verbose_info.figure_handles = [verbose_info.figure_handles, debug_image_handle];
end;

return;


%% ------------------------------------------------------------------------
function refined_corner = ICG_imposeLocalPlanarity (image, initial_point, neighbor_points, verbose)


current_corner_guess = [initial_point, neighbor_points];
current_corner_guess(3,:) = 1;
corners_pix = [0.5  0.5 0.5 0   1;...
               0.5  0   1   0.5 0.5];
corners_pix(3,:) = 1;

initial_homography = ICG_ransacFitHomography(corners_pix, current_corner_guess, 1e-1);

corner_point = ICG_normalizePoints(initial_homography * corners_pix(:,1));
refined_corner = corner_point(1:2);


%% --------------------------------------------------------------------------
function refined_corner = ICG_refineCheckerFrameProjective (image, initial_point, neighbor_points, template_type, verbose, USE_IPP)

num_template_cols_eighth = 12;
num_template_rows_eighth = 12;
num_template_rows = num_template_rows_eighth*8;
num_template_cols = num_template_cols_eighth*8;

persistent template t_type t_nre t_nce;
if ~isempty(template) && strcmp(template_type, t_type) && t_nre == num_template_rows_eighth && t_nce == num_template_cols_eighth,
else,
    template = ICG_createRegistrationTemplate(template_type, num_template_rows_eighth*8, num_template_cols_eighth*8, 0);
    t_type = template_type;
    t_nre = num_template_rows_eighth;
    t_nce = num_template_cols_eighth;
end;

current_corner_guess = [initial_point, neighbor_points];
current_corner_guess(3,:) = 1;

try,
    initial_homography = ICG_fitAffinity2D(template.neighbors, current_corner_guess);
catch,
    refined_corner = [];
    return;
end;

% Check whether region of interest is completely inside image
safety_delta = 20;
working_region = [-safety_delta                  -safety_delta;
                  num_template_cols+safety_delta -safety_delta;
                  -safety_delta                  num_template_rows+safety_delta;
                  num_template_cols+safety_delta num_template_rows+safety_delta]';
working_region(3,:) = 1;
transformed_working_region = round(ICG_normalizePoints(initial_homography * working_region));
transformed_midpoints = round(ICG_normalizePoints(initial_homography * template.feature_points(:,1:2)));

[num_rows, num_cols] = size(image);
inside_mask = ICG_pointsInsideRegion(transformed_working_region(1:2,:), [1 num_cols; 1 num_rows]);
if any(~inside_mask),
    refined_corner = [];
    return;
end;

sub_image_intervals = ICG_boundingBox(transformed_working_region(1:2,:));
sub_image = image(sub_image_intervals(2,1):sub_image_intervals(2,2), sub_image_intervals(1,1):sub_image_intervals(1,2));
t_img_to_sub = eye(3);
t_img_to_sub(1:2,3) =  [1-sub_image_intervals(1,1); 1-sub_image_intervals(2,1)];
initial_homography = t_img_to_sub * initial_homography;

initial_homography = initial_homography / initial_homography(3,3);
initial_homography(1,1) = initial_homography(1,1)-1;
initial_homography(2,2) = initial_homography(2,2)-1;

save_pattern = template.pattern;
is_inverted = false;
if strcmp(template_type, 'checkerboard'),
    template.pattern = invertTemplateIfNecessary(template.pattern, image, transformed_midpoints);
    is_inverted = true;
end;

filter = fspecial('gaussian', 5, 1.5);
% sub_image = imfilter(sub_image, filter, 'same', 'replicate');
sub_image = sub_image / max(sub_image(:));
pattern = imfilter(double(template.pattern), filter, 'same', 'replicate');

maxiter = 300;
num_iter = 1;
while 1,
    if USE_IPP,
        fit = ICG_projectiveRegistrationIPP (sub_image*255, pattern*255, initial_homography, 10, 0);
    else
        fit = ICG_projectiveRegistration (sub_image*255, pattern*255, initial_homography, 10, 0);
    end;
    [min_error, min_index] = min(cat(1,fit.rms_error));
    
    if min_index >= 10 && num_iter <= maxiter,
        initial_homography = fit(min_index(1)).warp_p;
    else
        refined_homography = fit(min_index(1)).warp_p;
        break;
    end;
    num_iter = num_iter+1;
end;

if verbose >= 3,
    fprintf(1, '\tMin Error: %.2f, Num Iterations: %.d\n', min_error, num_iter*10);
end;

refined_homography(1,1) = refined_homography(1,1) +1;
refined_homography(2,2) = refined_homography(2,2) +1;
refined_homography = inv(t_img_to_sub) * refined_homography;

corner_point = ICG_normalizePoints(refined_homography * template.center_xy);
refined_corner = corner_point(1:2);

if is_inverted
    template.pattern = save_pattern;
end


%% --------------------------------------------------------------------------
function refined_corner = ICG_refineCheckerFrameSaddle (image, initial_point, neighbor_points, template_type, verbose, USE_IPP)

num_template_cols_eighth = 12;
num_template_rows_eighth = 12;
num_template_rows = num_template_rows_eighth*8;
num_template_cols = num_template_cols_eighth*8;

persistent template t_type t_nre t_nce;
if ~isempty(template) && strcmp(template_type, t_type) && t_nre == num_template_rows_eighth && t_nce == num_template_cols_eighth,
else
    template = ICG_createRegistrationTemplate(template_type, num_template_rows_eighth*8, num_template_cols_eighth*8, 0);
    t_type = template_type;
    t_nre = num_template_rows_eighth;
    t_nce = num_template_cols_eighth;
end;

current_corner_guess = [initial_point, neighbor_points];
current_corner_guess(3,:) = 1;

try
    initial_homography = ICG_fitAffinity2D(template.neighbors, current_corner_guess);
catch
    refined_corner = [];
    return;
end;

% Check whether region of interest is completely inside image
safety_delta = 20;
working_region = [-safety_delta                  -safety_delta;
                  num_template_cols+safety_delta -safety_delta;
                  -safety_delta                  num_template_rows+safety_delta;
                  num_template_cols+safety_delta num_template_rows+safety_delta]';
working_region(3,:) = 1;
transformed_working_region = round(ICG_normalizePoints(initial_homography * working_region));
transformed_midpoints = round(ICG_normalizePoints(initial_homography * template.feature_points(:,1:2)));

[num_rows, num_cols] = size(image);
inside_mask = ICG_pointsInsideRegion(transformed_working_region(1:2,:), [1 num_cols; 1 num_rows]);
if any(~inside_mask),
    refined_corner = [];
    return;
end;

sub_image_intervals = ICG_boundingBox(transformed_working_region(1:2,:));
sub_image = image(sub_image_intervals(2,1):sub_image_intervals(2,2), sub_image_intervals(1,1):sub_image_intervals(1,2));
t_img_to_sub = eye(3);
t_img_to_sub(1:2,3) =  [1-sub_image_intervals(1,1); 1-sub_image_intervals(2,1)];
initial_homography = t_img_to_sub * initial_homography;

initial_homography = initial_homography / initial_homography(3,3);

if strcmp(template_type, 'checkerboard'),
    template.pattern = invertTemplateIfNecessary(template.pattern, image, transformed_midpoints);
end;

filter = fspecial('gaussian', 5, 1.5);
% sub_image = imfilter(sub_image, filter, 'same', 'replicate');
sub_image = sub_image / max(sub_image(:));

warped = ICG_projectiveWarpIPP(im2uint8(sub_image), inv(initial_homography), [1, 1, size(sub_image, 2), size(sub_image, 1)], [1, 1, size(template.pattern, 2), size(template.pattern,1)]);
warped_f = imfilter(warped, fspecial('gaussian', 15, 9), 'same', 'replicate');

saddle_region = warped_f(31:end-30, 31:end-30);
[xg, yg] = meshgrid(31:size(warped_f,2)-30, 31:size(warped_f,1)-30);

[r, c] = size(saddle_region);
cx = mean(xg(:));
cy = mean(yg(:));
ptn = [xg(:)-cx, yg(:)-cy]';
scale = sqrt(2) / mean(sqrt(sum(ptn.^2, 1)));
ptn = ptn*scale;
A = [ptn(1,:)'.^2, ptn(1,:)'.*ptn(2,:)', ptn(2,:)'.^2, ptn(1,:)', ptn(2,:)', ones(size(ptn,2),1)];
b = double(saddle_region(:))/255;
surface = A\b;
saddle_n = -inv([[2*surface(1), surface(2)]; [surface(2) 2*surface(3)]]) * [surface(4); surface(5)];
saddle = saddle_n / scale + [cx; cy];
saddle(3) = 1;

refined_corner = inv(t_img_to_sub)*initial_homography*saddle;
refined_corner(3) = [];


%% --------------------------------------------------------------------------
function refined_corner = ICG_refineCheckerFrameCirclesCorr (image, initial_point, neighbor_points, template_type, verbose, USE_IPP)

num_template_cols_eighth = 12;
num_template_rows_eighth = 12;
num_template_rows = num_template_rows_eighth*8;
num_template_cols = num_template_cols_eighth*8;

filter0 = fspecial('gaussian', 3, 0.5);
filter1 = fspecial('gaussian', 11, 3);
filter2 = fspecial('sobel');
filter3 = filter2';

persistent template t_type t_nre t_nce corr_template;
if ~isempty(template) && strcmp(template_type, t_type) && t_nre == num_template_rows_eighth && t_nce == num_template_cols_eighth,
else,
    template = ICG_createRegistrationTemplate(template_type, num_template_rows_eighth*8, num_template_cols_eighth*8, 0);
    t_type = template_type;
    t_nre = num_template_rows_eighth;
    t_nce = num_template_cols_eighth;
    corr_template = template.pattern(template.feature_points(2,1)-num_template_rows_eighth-8:template.feature_points(2,1)+num_template_rows_eighth+8, ...
                                     template.feature_points(1,1)-num_template_cols_eighth-8:template.feature_points(1,1)+num_template_cols_eighth+8);
    corr_template = imfilter(corr_template, filter1, 'same', 'replicate');
%     corr_template_h = imfilter(corr_template, filter2, 'same', 'replicate');
%     corr_template_v = imfilter(corr_template, filter3, 'same', 'replicate');
%     corr_template = corr_template_h.^2 + corr_template_v.^2;
    corr_template = corr_template(9:end-8, 9:end-8);
end;

current_corner_guess = [initial_point, neighbor_points];
current_corner_guess(3,:) = 1;

try,
    initial_homography = ICG_fitAffinity2D(template.neighbors, current_corner_guess);
catch,
    refined_corner = [];
    return;
end;

% Check whether region of interest is completely inside image
safety_delta = 20;
working_region = [-safety_delta                  -safety_delta;
                  num_template_cols+safety_delta -safety_delta;
                  -safety_delta                  num_template_rows+safety_delta;
                  num_template_cols+safety_delta num_template_rows+safety_delta]';
working_region(3,:) = 1;
transformed_working_region = round(ICG_normalizePoints(initial_homography * working_region));
transformed_midpoints = round(ICG_normalizePoints(initial_homography * template.feature_points(:,1:2)));

[num_rows, num_cols] = size(image);
inside_mask = ICG_pointsInsideRegion(transformed_working_region(1:2,:), [1 num_cols; 1 num_rows]);
if any(~inside_mask),
    refined_corner = [];
    return;
end;

sub_image_intervals = ICG_boundingBox(transformed_working_region(1:2,:));
sub_image = image(sub_image_intervals(2,1):sub_image_intervals(2,2), sub_image_intervals(1,1):sub_image_intervals(1,2));
t_img_to_sub = eye(3);
t_img_to_sub(1:2,3) =  [1-sub_image_intervals(1,1); 1-sub_image_intervals(2,1)];
initial_homography = t_img_to_sub * initial_homography;


refined_homography = initial_homography;
% corner_history = {};
for iter = 1:3,
    if USE_IPP,
        [wimg status] = ICG_projectiveWarpIPP(im2uint8(sub_image), inv(refined_homography), ...
            [1 1 size(sub_image, 2) size(sub_image, 1)], ...
            [1 1 size(template.pattern, 2) size(template.pattern, 1)]);
    else,
        [wimg status] = ICG_projectiveWarpIPP(im2uint8(sub_image), inv(refined_homography), ...
            [1 1 size(sub_image, 2) size(sub_image, 1)], ...
            [1 1 size(template.pattern, 2) size(template.pattern, 1)]);
    end;

%     wimg = imfilter(wimg, filter0);
    wimg = single(wimg);
%     filtered_warped_h = imfilter(wimg, filter2);
%     filtered_warped_v = imfilter(wimg, filter3);
%     filtered_warped = filtered_warped_h.^2 + filtered_warped_v.^2;
%     filtered_warped = double(filtered_warped);
    filtered_warped = im2double(wimg);
    
%     filtered_warped(1:3,:) = 0;
%     filtered_warped(end-2:end,:) = 0;
%     filtered_warped(:,1:3) = 0;
%     filtered_warped(:,end-2:end) = 0;

    corr_pattern = ICG_normXCorr2 (corr_template, filtered_warped, 'same');
    [rows, cols] = size(corr_pattern);
    offset_xy = [cols/2+1, rows/2+1]-1;
    vxg = 1:cols/2;
    vyg = 1:rows/2;
    [xg, yg] = meshgrid(vxg, vyg);
    corners_subpix = ones(3, 4);
    for i=1:2,
        for j=1:2,
            quadrant = corr_pattern(vyg+(i-1)*offset_xy(2), vxg+(j-1)*offset_xy(1));
            [maxval, maxind] = max(quadrant(:));
            temp_xg = xg+(j-1)*offset_xy(1);
            temp_yg = yg+(i-1)*offset_xy(2);
            max_x = temp_xg(maxind);
            max_y = temp_yg(maxind);

            [x_offset, y_offset, max_f] = ICG_findPeak(quadrant,1);
            corners_subpix(1:2, (i-1)*2+(j-1)+1) = [max_x+x_offset; max_y+y_offset];
        end; 
    end;
    
    
    
    corners_subpix = corners_subpix(:, [1 2 4 3]);
    corners_subpix = ICG_normalizePoints(refined_homography * corners_subpix);
%     corner_history{end+1} = corners_subpix;
    refined_homography = ICG_fitHomography2D(template.feature_points, corners_subpix);
end;

refined_homography = inv(t_img_to_sub) * refined_homography;

corner_point = ICG_normalizePoints(refined_homography * template.center_xy);
refined_corner = corner_point(1:2);


%% --------------------------------------------------------------------------
function refined_corner = ICG_refineCheckerFrameCircles (initial_point, neighbor_points, circle_centers, template, verbose);

% For each point, estimate homography, warp circles, find nearest segmented
% circles, re-estimate homography on circles, map center point.

current_corner_guess = [initial_point, neighbor_points];
current_corner_guess(3,:) = 1;
try
    initial_homography = ICG_fitAffinity2D(template.neighbors, current_corner_guess);
catch
    refined_corner = [];
    return;
end;

circle_centers_template = template.feature_points;
% circle_centers_template(3,:) = 1;
transformed_circle_centers = ICG_normalizePoints(initial_homography * circle_centers_template);

% Find circle centers closest to transformed guesses
refined_center = zeros(3, size(transformed_circle_centers, 2))*NaN;
for i=1:size(transformed_circle_centers, 2),
    current_center = transformed_circle_centers(:, i);
    sqr_distances = sum((circle_centers(1,:)-current_center(1)).^2 + (circle_centers(2,:)-current_center(2)).^2, 1);
    [min_dist, min_ind] = min(sqr_distances);
    if ~isempty(min_ind) && sqrt(min_dist) < 3,
        refined_center(:, i) = circle_centers(:, min_ind);
    end;
end;
bad = isnan(refined_center(1,:));
if isempty(bad) || any(bad),
    refined_corner = [];
    return;
end;

center_xy = template.center_xy;
% center_xy(3) = 1;
refined_homography = ICG_fitHomography2D(circle_centers_template, refined_center);
refined_corner = ICG_normalizePoints(refined_homography * center_xy);
refined_corner(3) = [];


%% ------------------------------------------------------------------------
function refined_corner = ICG_refineCheckerFrameLinecross (image, initial_point, neighbor_points, verbose);
% refined_corner = ICG_refineCheckerFrameLinecross (image, initial_point, neighbor_points, verbose);

sigma = 2;

up_vector = neighbor_points(:,1) - neighbor_points(:,2);
right_vector = neighbor_points(:,4) - neighbor_points(:,3);
frame_orientation_up = atan2(up_vector(2), up_vector(1));
frame_orientation_right   = atan2(right_vector(2), right_vector(1));

up_norm = norm(up_vector);
right_norm = norm(right_vector);
if up_norm<sigma*6 || right_norm<sigma*6,
    refined_corner = [];
    return;
end;

initial_up(:,1) = neighbor_points(:,2) + up_vector / 6;
[up_edge_pos_1, up_edge_strength_1] = ICG_profileEdge (image, initial_up(:,1), frame_orientation_right, right_norm/4, sigma);

initial_up(:,2) = neighbor_points(:,2) + up_vector / 3;
[up_edge_pos_2, up_edge_strength_2] = ICG_profileEdge (image, initial_up(:,2), frame_orientation_right, right_norm/4, sigma);

initial_up(:,3) = neighbor_points(:,2) + up_vector / 3 * 2;
[up_edge_pos_3, up_edge_strength_3] = ICG_profileEdge (image, initial_up(:,3), frame_orientation_right, right_norm/4, sigma);

initial_up(:,4) = neighbor_points(:,2) + up_vector / 6 * 5;
[up_edge_pos_4, up_edge_strength_4] = ICG_profileEdge (image, initial_up(:,4), frame_orientation_right, right_norm/4, sigma);


initial_right(:,1) = neighbor_points(:,3) + right_vector / 6;
[right_edge_pos_1, right_edge_strength_1] = ICG_profileEdge (image, initial_right(:,1), frame_orientation_up, up_norm/4, sigma);

initial_right(:,2) = neighbor_points(:,3) + right_vector / 3;
[right_edge_pos_2, right_edge_strength_2] = ICG_profileEdge (image, initial_right(:,2), frame_orientation_up, up_norm/4, sigma);

initial_right(:,3) = neighbor_points(:,3) + right_vector / 3 * 2;
[right_edge_pos_3, right_edge_strength_3] = ICG_profileEdge (image, initial_right(:,3), frame_orientation_up, up_norm/4, sigma);

initial_right(:,4) = neighbor_points(:,3) + right_vector / 6 * 5;
[right_edge_pos_4, right_edge_strength_4] = ICG_profileEdge (image, initial_right(:,4), frame_orientation_up, up_norm/4, sigma);

upper_edge_points = up_edge_pos_1;
upper_edge_points = [upper_edge_points, up_edge_pos_2];
upper_edge_points = [upper_edge_points, up_edge_pos_3];
upper_edge_points = [upper_edge_points, up_edge_pos_4];

right_edge_points = right_edge_pos_1;
right_edge_points = [right_edge_points, right_edge_pos_2];
right_edge_points = [right_edge_points, right_edge_pos_3];
right_edge_points = [right_edge_points, right_edge_pos_4];

if size(upper_edge_points, 2) < 4 || size(right_edge_points, 2) < 4,
    refined_corner = [];
    return;
end;


% if 1,
%     hold on, plot ([upper_points(1,1) upper_edge_pos(1,1)], [upper_points(2,1) upper_edge_pos(2,1)], '-b.');
%     hold on, plot ([upper_points(1,2) upper_edge_pos(1,2)], [upper_points(2,2) upper_edge_pos(2,2)], '-b.');
% 
%     hold on, plot ([lower_points(1,1) lower_edge_pos(1,1)], [lower_points(2,1) lower_edge_pos(2,1)], '-b.');
%     hold on, plot ([lower_points(1,2) lower_edge_pos(1,2)], [lower_points(2,2) lower_edge_pos(2,2)], '-b.');
% 
%     hold on, plot ([left_points(1,1) left_edge_pos(1,1)], [left_points(2,1) left_edge_pos(2,1)], '-b.');
%     hold on, plot ([left_points(1,2) left_edge_pos(1,2)], [left_points(2,2) left_edge_pos(2,2)], '-b.');
% 
%     hold on, plot ([right_points(1,1) right_edge_pos(1,1)], [right_points(2,1) right_edge_pos(2,1)], '-b.');
%     hold on, plot ([right_points(1,2) right_edge_pos(1,2)], [right_points(2,2) right_edge_pos(2,2)], '-b.');
% end;


[vertical_line, consensus, inliers, vertical_error] = ICG_ransacFitLine(upper_edge_points, 1);
if size(inliers, 2) < 3,
    refined_corner = [];
    return;
end;
   
[horizontal_line, consensus, inliers, horizontal_error] = ICG_ransacFitLine(right_edge_points, 1);
if size(inliers, 2) < 3,
    refined_corner = [];
    return;
end;

error = [vertical_error; horizontal_error];
refined_corner = ICG_normalizePoints(cross(vertical_line, horizontal_line));
refined_corner(3) = [];

%% ------------------------------------------------------------------------
function template = invertTemplateIfNecessary(template, image, transformed_midpoints);

if image(transformed_midpoints(2,1), transformed_midpoints(1,1)) > image(transformed_midpoints(2,2), transformed_midpoints(1,2)),
    template = ~template;
end;
