SPT
single particle tracking routines
==

The duration of this study spanned from Sep 2013 to Sep 2014, by Tzu-Chi Yen, under the supervision of Dr. Chia-Lung Hsieh.

The package consists of onw part:
    (I)	    [code-tracking]: General codes for SPT
 
#=============================================================
(I) [code-tracking]:
(I-i) tracking_main: Starting with main code *run.m* and adapting changes in *track.m* and *sif2tif.m*, one will be able to acquire trajectory/msd info from a series of properly named *.sif files.

(I-ii) localization_benchmark: Using *example_run.m*, you will obtain simulated images and their localization results. *make_multiSNr_plots.m* will help to plot the errors w.r.t. to SNr’s. Refer *Parthasarathy_readme.pdf* for more detailed explanation. 

(I-iii) OSS_benchmark: *crop_img.m* is main code for generating the time-resolved intensity fluctuation profile in a series of OSS benchmark experiment.

(I-iv) Post_processing:
	(I-iv-1) brownian_comparison: statistical analysis with simulated Brownian motion
	(I-iv-2) moving_D_analysis: Analysis for transient diffusion constant
	(I-iv-3) non-gaussian evolution analysis: Analysis for non-Gaussianity 
	(I-iv-4) velocity_autocorrelation: Trajectory velocity autocorrelation
	(I-iv-5) miscellaneous: Hurst exponent, waiting-time distribution, etc
