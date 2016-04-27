--------------------------------------------------------------------------------
| Automatic feature extraction for camera calibration                          |
| D. Ferstl, C. Reinbacher, G. Riegler, M. Rüther, H. Bischof                  |
| Graz University of Technology, Institute for Computer Vision and Graphics    |
| (c) 2015                                                                     |
--------------------------------------------------------------------------------

This software accompanies our paper

@inproceedings{ferstl2015,
  author = {David Ferstl and Christian Reinbacher and Gernot Riegler and 
            Matthias Ruether and Horst Bischof},
  title = {Learning Depth Calibration of Time-of-Flight Cameras},
  booktitle = {Proceedings of British Machine Vision Conference, (BMVC)},
  year = {2015},
  month = {September},
} 

The presented calibration target and automatic feature extraction are not
limited to depth cameras but can also be used for conventional cameras. The
provided code is designed to be used as an addon to the widely known camera
calibration toolbox of [1].

The used calibration target consists of a central marker and circular feature
points. Our automatic feature detection starts by searching for the central
marker and then iteratively refining the circular markers around the central
marker. Compared to standard checkerboard targets, our methods has the following
advantages:
- Target does not have to be visible as a whole
- Detection of groups of circular patterns is more robust to perspective 
  distortions than line crossings
- Feature detection is more accurate for low-resolution cameras (like ToF or 
  Event Cameras)

--------------------------------------------------------------------------------
| How to use the code                                                          |
--------------------------------------------------------------------------------
The provided source code should be used as an addon the the Bouguet Camera 
Calibration Toolbox. Installation therefore amounts to:
- Downloading of the latest toolbox from [1]
- Downloading and extracting of our addon in the same folder
- Running autocalibration.m and selecting the images from testdata/image_xxx.jpg 
  starts the mono calibration of the camera.
- The calibration target can be created using the make_target.m function. 
  Remember to measure it after printing!
- To use it with the GUI of the Toolbox, simply start calib_gui_normal_auto.m 
  which asks for the target parameters interactively.
- Stereo calibration requires the use of calib_stereo_auto.m instead of 
  calib_stereo.m because our method does not detect all grid points in all 
  images!

The following parameters have to be set:
- parameters.approx_marker_width_pixels: Approximate minimum size of the center 
  marker in pixels
- parameters.grid_width_mm: Grid width (distance between points) in millimeters
- parameters.checker_aspect_ratio: Aspect ratio (= height/width)
- parameters.grid_coordinates_h:Horizontal grid dimensions (i.e. -11:11)
- parameters.grid_coordinates_v: Vertical grid dimensions (i.e. -18:16)

The target can be created by using the function function template = make_target
(grid_width_pixels, grid_width_mm, grid_coordinates_h, grid_coordinates_v),
i.e.: target = make_target(240,5,-18:18,-10:10);

--------------------------------------------------------------------------------
| Version History                                                              |
--------------------------------------------------------------------------------
- v 0.1: Initial Release


--------------------------------------------------------------------------------
| References                                                                   |
--------------------------------------------------------------------------------

[1]: http://www.vision.caltech.edu/bouguetj/calib_doc/

The inverse compositional algorithm has been taken from:
--------------------------------------------------------------------------------
URL: http://www.ri.cmu.edu/projects/project_515.html

Bibtex reference:

@article{Baker_2004_4293,
	author = "Simon Baker and Iain Matthews",
	title = "Lucas-Kanade 20 Years On: A Unifying Framework Part 1: The Quantity 
	Approximated, the Warp Update Rule, and the Gradient Descent Approximation",
	journal = "International Journal of Computer Vision",
	year = "2004"
}

The used files are:
- hessian.m
- homo_ic.m
- init_h.m
- jacobian_h.m
- quadtobox_h.m
- sd_images.m
- sd_update.m
- warp_h.m

The feature detection is based on ARToolkitPlus, taken from:
--------------------------------------------------------------------------------
URL: http://studierstube.icg.tugraz.at/handheld_ar/artoolkitplus.php

@inproceedings{Wagner_2007,
	author = {Daniel Wagner and Dieter Schmalstieg},
	title =  {ARToolKitPlus for Pose Tracking on Mobile Devices},
	booktitle = {Proceedings of 12th Computer Vision Winter Workshop (CVWW'07)}
	month = {February},
	year = {2007},
}

The normalized cross correlation from OpenCV is used:
--------------------------------------------------------------------------------

The source code for the wrapper is in the normxcorr2_src/ subdirectory and 
should compile out of the box when using CMake, OpenCV and Matlab under both 
windows and linux.

For those who do not want to use CMake, the following MATLAB command should do:
(example for Visual Studio 2013 on a x64 machine)
 mex -I'path\to\opencv\build\include\opencv' ...
     -I'path\to\opencv\build\include' ...
     path\to\opencv\build\x64\vc12\lib\opencv_core2411.lib ...   
     path\to\opencv\build\x64\vc12\lib\opencv_imgproc2411.lib ...
     normxcorr2_mex.cpp
