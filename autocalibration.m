% Autocalibration addon for Bouguet Toolbox version July 2014

toolbox_path = fullfile(pwd,'TOOLBOX_calib/');

if(~exist(toolbox_path, 'dir'))
    !wget http://www.vision.caltech.edu/bouguetj/calib_doc/download/toolbox_calib.zip
    !unzip toolbox_calib.zip
    !rm toolbox_calib.zip
end

addpath(toolbox_path);

image_dir = './testdata';
%% Read in images
path__ = pwd;
cd(image_dir);
data_calib;
ima_read_calib;
cd(path__);clear path__
%% Set parameters
parameters.approx_marker_width_pixels = 60;
parameters.grid_width_mm = 10.218;
parameters.image_threshold = [];
parameters.checker_aspect_ratio = 10.18/10.218;
parameters.corner_extraction_method = 'projective';
parameters.target_type = 'circles';
parameters.template_type = 'circles';
parameters.grid_coordinates_h = -11:11;
parameters.grid_coordinates_v = -18:16;
parameters.pixel_size_wh_mm = [6, 6]*1e-3;
parameters.feature_subsampling_rate = 2;
parameters.USE_IPP=0;
parameters.verbose=1;

%% Run feature extraction
auto_calib;

%% Run calibration
go_calib_optim;