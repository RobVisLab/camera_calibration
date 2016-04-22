% Automatic feature extraction
[folder, name, ext] = fileparts(mfilename('fullpath'));
addpath(fullfile(folder,'ICG_functions/'));

if(~exist('parameters','var'))
    % Set default parameters
    parameters.image_threshold = [];
    parameters.corner_extraction_method = 'projective';
    parameters.target_type = 'circles';
    parameters.template_type = 'circles';
    parameters.pixel_size_wh_mm = [6, 6]*1e-3;
    parameters.feature_subsampling_rate = 2;
    parameters.USE_IPP=0;
    parameters.verbose=1;
    
    % Ask for mandatory parameters
    parameters.approx_marker_width_pixels = input('Enter approx. minimal marker width in pixels: ');
    parameters.grid_width_mm = input('Enter grid with in millimeters: ');
    parameters.checker_aspect_ratio = input('Enter grid aspect ratio (= height/width): ');
    parameters.grid_coordinates_h = input('Enter horizontal grid dimensions (i.e. -11:11): ');
    parameters.grid_coordinates_v = input('Enter vertical grid dimensions (i.e. -18:16): ');
end
dX = parameters.grid_width_mm;
dY = parameters.grid_width_mm*parameters.checker_aspect_ratio;
verbose_info.verbose = parameters.verbose;
verbose_info.figure_handles = [];
template = ICG_createCalibrationMarker(16*7, 7, parameters.template_type, 0);
template_pattern = template.pattern(template.corner_points(2,1):template.corner_points(2,3), template.corner_points(1,1):template.corner_points(1,3));
template_pattern = imresize(template_pattern, [32 32], 'bilinear');

corner_points = cell(1, n_ima);%{};
grid_coords = cell(1, n_ima);%{};
image_stack = cell(1, n_ima);
for i=1:n_ima,
    image_stack{i} = eval(['I_' num2str(i)]);
end

for i=1:n_ima,
    fprintf('Processing image %d/%d', i, n_ima);
    localverbose_info = verbose_info;
    [active_images(i), corner_points{i}, grid_coords{i}] = ICG_extractCornersLoop (image_stack{i}, template_pattern, parameters, localverbose_info,i);
end

for i=1:n_ima,
    if active_images(i)
        eval(['X_' num2str(i) '=grid_coords{i};']);
        eval(['x_' num2str(i) '=corner_points{i}(1:2,:);']);
    else
        eval(['X_' num2str(i) '=nan(3,1);']);
        eval(['x_' num2str(i) '=nan(2,1);']);
    end
end
