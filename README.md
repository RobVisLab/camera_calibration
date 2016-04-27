<html>
<head>
	<meta http-equiv="content-type" content="text/html; charset=utf-8" />
	<link href="readme/style.css" rel="stylesheet" type="text/css" media="screen" />
</head> 
<body>
	<hr />
	<div id="page">
		<div id="content">
			<div class="post">
				<h1 class="title">Automatic Camera Calibration</h1>
				<div class="entry">
					This page accompanies our paper <a target="_blank" href="https://rvlab.icg.tugraz.at/documents/ferstl/bmvc15_final.pdf">[1]</a> on
					automatic calibration of depth cameras. The presented calibration
					target and automatic feature extraction are not limited to depth
					cameras but can also be used for conventional cameras. The provided
					code is designed to be used as an addon to the widely known camera
					calibration toolbox of <a
					href="http://www.vision.caltech.edu/bouguetj/calib_doc/index.html">Jean-Yves
					Bouguet</a>
				</div>
				<h2 class="title">Calibration Target</h2>
				<div class="entry">
					The used calibration target consists of a central marker and
					circular patterns:<br></br> <img alt="Calibration Target"
					src="https://rvlab.icg.tugraz.at/project_page/project_tofusion/calibration/struc_patt.png" width=300px><br></br> Our
					automatic feature detection starts by searching for the central
					marker and then iteratively refining the circular markers around
					the central marker (depticted as black dashed line). Compared to
					standard checkerboard targets, our methods has the following
					advantages:
					<ul>
						<li>Target does not have to be visible as a whole</li>
						<li>Detection of groups of circular patterns is more robust
							to perspective distortions than line crossings</li>
						<li>Feature detection is more accurate for low-resolution
							cameras (like ToF or Event Cameras)</li>
					</ul>
				</div>
				<h2 class="title">Example Detection Result</h2>
				<div class="entry">
					The following images show the result of the automatic feature
					detection on two exemplary images from the paper. The calibration
					target is detected in the gray value image and reprojected to the
					corresponding depth image from a Microsoft Kinect v2.0: <br></br> <img
						alt="Detection result gray image"	src="https://rvlab.icg.tugraz.at/project_page/project_tofusion/calibration/detection_rgb.jpg" height=270px> <img
						alt="Detection result projected to depth image"	src="https://rvlab.icg.tugraz.at/project_page/project_tofusion/calibration/projection_depth.jpg" height=270px><br></br>
				</div>
			</div>
			<br><br>
			<div class="post">
			<h2 class="title">How to use the code</h2>
				<div class="entry">
					The provided source code is used as an addon the the Bouguet
					Camera Calibration Toolbox. Installation therefore amounts to:
					<ul>
						<li>Clone the calibration toolbox from GitHub <a
							href="https://github.com/RobVisLab/camera_calibration">[GitLab Link]</a></li>
						<li>With the <i>download_script.sh</i> the Bouguet Toolbox can be updated to the latest version</li>
						<li>Running <i>autocalibration.m</i> and selecting the images
							from <i>testdata/image_xxx.jpg</i> starts the mono calibration of
							the camera.
						</li>
						<li>The calibration target can be created using the <i>make_target.m</i> function. Remember to measure it after printing! </li>
						<li>To use it with the GUI of the Toolbox, simply start <i>calib_gui_normal_auto.m</i>
							which asks for the target parameters interactively.</li>
						<li>Stereo calibration requires the use of <i>calib_stereo_auto.m</i> instead of <i>calib_stereo.m</i> because our method does not detect all grid points in all images! </li>
					</ul>
					The following parameters have to be set:
					<ul>
						<li><i>parameters.approx_marker_width_pixels:</i> Approximate
							minimum size of the center marker in pixels</li>
						<li><i>parameters.grid_width_mm:</i> Grid width (distance
							between points) in millimeters</li>
						<li><i>parameters.checker_aspect_ratio:</i> Aspect ratio (=
							height/width)</li>
						<li><i>parameters.grid_coordinates_h:</i>Horizontal grid
							dimensions (i.e. -11:11)</li>
						<li><i>parameters.grid_coordinates_v:</i> Vertical grid
							dimensions (i.e. -18:16)</li>
					</ul>
					The target can be created by using the function <i>function template = make_target (grid_width_pixels, grid_width_mm, grid_coordinates_h, grid_coordinates_v)</i>, i.e.:
					<i>target = make_target(240,5,-18:18,-10:10);</i>
				</div>
			</div>
			</div>
			<br><br>
			<div class="post">
				<h2 class="title">How to cite the materials</h2>
				<div class="entry">
					We grant permission to use the code on this website. If you if you
					use the code in your own publication, we request that you cite our
					paper <a href="#publications">[1]</a>. If you want to cite this
					website, please use the URL
					"https://rvlab.icg.tugraz.at/calibration/".
				</div>
			</div>
			<br><br>
			<div class="publications" id="publications">
				<h2 class="title">References</h2>
					<ol>
					<li>
					 <a target="_blank" target="_blank" class="title" href="https://rvlab.icg.tugraz.at/documents/ferstl/bmvc15_final.pdf">Learning Depth Calibration of Time-of-Flight Cameras</a>
					&nbsp;
					  <a target="video" href="https://youtu.be/RM78k8M2qiw" ><img src="https://rvlab.icg.tugraz.at/images/youtube.png" height=15px/></a>
					  <a target="_blank" class="title" href="https://rvlab.icg.tugraz.at/documents/ferstl/bmvc15_poster.pdf">[supp]</a> <br>
					<span class="authors"><a href=https://rvlab.icg.tugraz.at/personal_page/personal_page_ferstl.html>David Ferstl</a>, <a href=https://rvlab.icg.tugraz.at/personal_page/personal_page_christian.htm>Christian Reinbacher</a>, <a href=https://rvlab.icg.tugraz.at/personal_page/personal_page_gernot.html>Gernot Riegler</a>, Matthias Ruether, and <a href=http://www.icg.tugraz.at/Members/author/bischof target=_blank>Horst Bischof</a></span><br>
					In <span class="in">Proceedings of British Machine Vision Conference, (BMVC)</span>, <span class="year">2015</span> <br/> Abstract </a><br />
					<div class="abstract" id="ferstl2015a" style="display: none">
					We present a novel method for an automatic calibration of modern consumer
					 Time-of-Flight cameras. Usually, these sensors come equipped with an integrated color
					 camera. Albeit they deliver acquisitions at high frame rates they usually suffer
					 from incorrect calibration and low accuracy due to multiple error sources. Using
					 information from both cameras together with a simple planar target, we will show
					 how to accurately calibrate both color and depth camera and tackle most error
					 sources inherent to Time-of-Flight technology in a unified calibration framework.
					 Automatic feature detection minimizes user interaction during calibration. We
					 utilize a Random Regression Forest to optimize the manufacturer supplied depth 
					 measurements. We show the improvements to commonly used depth calibration
					 methods in a qualitative and quantitative evaluation on multiple scenes acquired
					 by an accurate reference system for the application of dense 3D reconstruction.</div></li>
					</ol>
			</div>
		</div>
	</div>
</body>
</html>
