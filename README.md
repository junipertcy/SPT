single particle tracking routines
==================

The duration of this study spanned from Sep 2013 to Sep 2014, by Tzu-Chi Yen, under the supervision of Dr. Chia-Lung Hsieh.

The package components
==================
* [code-tracking]: General codes for SPT
	* tracking_main: Starting with main code run.m and adapting changes in track.m and sif2tif.m, one will be able to acquire trajectory/msd info from a series of properly named *.sif files.

	* localization_benchmark: Using example_run.m, you will obtain simulated images and their localization results. make_multiSNr_plots.m will help to plot the errors w.r.t. to SNr’s. Refer Parthasarathy_readme.pdf for more detailed explanation. 

	* OSS_benchmark: crop_img.m is main code for generating the time-resolved intensity fluctuation profile in a series of OSS benchmark experiment.

	* Post_processing:
		* brownian_comparison: statistical analysis with simulated Brownian motion
		* moving_D_analysis: Analysis for transient diffusion constant
		* non-gaussian evolution analysis: Analysis for non-Gaussianity 
		* velocity_autocorrelation: Trajectory velocity autocorrelation
		* miscellaneous: Hurst exponent, waiting-time distribution, etc
