% run_short.m
%
% This MATLAB script is for tracking of the immobile dyes for the benchmark
% of the oxygen scavenging system. Note in order to make particle detection
% sensitive, one need to set high enough *pth_int*, *L*, *t_cutoff* and
% *pth*, but make *k_rank* as close to the experiment as possible.
%
%
%
% Function calls:
%       tracker.m
%       detect_particles_IAMS_*.m
%       radialcenter.m
%       link_trajectories_IAMS.m   
%       sort_trajectories_IAMS.m
%       MSD_analyzer.m
% 
% Tzu-Chi Yen
% Feb. 17, 2014
%
% Project began Sep. 2013
%


tic
% ===== Begin User-defined Parameters =====
%str_path = '/Users/junipe/mnt/test_images/';
w = 5;        	 % radius of neighborhood: > particle radius, 
              	 % < interparticle distance
ww = 5;          % crop size of the image
pth = 0.1;       % upper intensity percentile for detecting
                 % particles (rough estimation)
pth_int = 580;   % least intensity bound for detecting
                 % particles (important as to find true dyes)

t_cutoff = 60; % cutoff length for truncation of trajectories

L = 2;            % maximum displacement between frames
k_rank = 10;      % maximum number of particles discriminated

N_bit = 16; 
pixel_size = 96e-3; % 48nm, v711, 20um pixel size, 75cn lens, 100x OBJ

img=[];
img=I_g(:,:,:,1);
    
%=============================================================
% Begin SPT calculation
%=============================================================
% First and foremost, the videos will be loaded and fed into a
% localization algorithm to extract the respective trajectories
    
    N_start = 1; % Start frame for localization
    N_end = 2000; % End frame for localization
    N_img = N_end-N_start+1;
    frame_rate = 100;
    frame_time = 1/frame_rate;
    tracker;   % save --append *peaks*
    MSD_analyzer;    % save --append *traj*
 
toc
