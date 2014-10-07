% run.m
%
% This MATLAB script is the main code for loading/processing a
% series of videos, and localizes their trajectories. 
%
%
% Function calls:
%       tracker.m
%       detect_particles_IAMS_*.m
%       radialcenter.m
%       link_trajectories_IAMS.m   
%       sort_trajectories_IAMS.m
%       MSD_calculator.m
%       MSD_analyzer.m
%       CDF_analyzer.m
% 
% Tzu-Chi Yen
% Feb. 17, 2014
%
% Project began Sep. 2013



tic
% ===== Begin User-defined Parameters =====
%str_path = '/Users/junipe/mnt/test_images/';
str_path = '/Users/junipe/Desktop/';

% Do make sure you have rename the video series into the format
% name_$number.sif or name_$number.tif
str_filename_main = 'PCA_PCD_Trolox_';
n_videos = 1;      % number of videos analysed

w = 5;        	 % radius of neighborhood: > particle radius, 
               	 % < interparticle distance
ww = 5;          % crop size of the image       
pth = 0.2;       % upper intensity percentile for detecting
                 % particles (rough estimation)
pth_int = 600;   % upper intensity *absolute* for detecting
                 % particles
L = 1;           % maximum displacement between frames
k_rank = 20;    % maximum number of particles discriminated
T_step = [1:20];

N_bit = 16; 
pixel_size = 96e-3; % 96 nm

%=============================================================
% parameters for saving the output
%=============================================================

str_method = 'r';
str_w = num2str(w);
str_pth = num2str(pth);
str_L = int2str(L);
str_k_rank = int2str(k_rank);
str_n_videos = int2str(n_videos);

savefile = ['/Users/junipe/Desktop/',str_method,'_SPT_',str_filename_main,'[',str_w,'_',str_k_rank,'_',str_pth,'_',str_L,'_',str_n_videos,']','.mat'];

save(savefile,'pth','w','ww','L','T_step','k_rank');

% ===== End User-defined Parameters =====

%=============================================================
% Begin SPT calculation
%=============================================================
% First and foremost, the videos will be loaded and fed into a
% localization algorithm to extract the respective trajectories


for i_spec=1:n_videos
    str_filename = sprintf('%s%d',str_filename_main,i_spec);

    % if *.sif file is loaded, you have to modify the codes in the
    % first section of tracker.m, ie function sif2tif should be called.
    filename = sprintf('%s%s%s',str_path,str_filename,'.sif');
    savevideo = [str_filename,'.tif'];
    
    
%    info = func_ext_uli_sifread2(filename);
%     Width = info.Width; 
%     Height = info.Height; 
%     resolution = info.resolution;
%     NumFrames = info.NumFrames;

    Width = 128;
    Height = 128;
    resolution = [128 128];
    NumFrames = 5000;
   
%     N_start = 1; % Start frame for localization
%     N_end = info.NumFrames; % End frame for localization
%     N_img = N_end-N_start+1;
%     frame_rate = info.framerate;
%     frame_time = 1/frame_rate;

    N_start = 1; % Start frame for localization
    N_end = 5000; % End frame for localization
    N_img = N_end-N_start+1;
    frame_rate = 100;
    frame_time = 1/frame_rate;
    tracker;   % save --append *peaks*
    MSD_analyzer;    % save --append *traj*

    
    
% clear(str_peaks,str_traj,'all_part','trajectories','peak','peaks','traj*')
end
load(savefile);

% Here, multi-video trajectories are recombined into a single variable
traj=[];
%traj_I=[];
for i=1:n_videos
traj=[traj eval(sprintf('%s%d','traj_',i))];
%traj_I=[traj_I eval(sprintf('%s%d','traj_I_',i))];
end
 
% % ===== BEGIN ensemble and time-averaged MSD calculation =====
[msd] = MSD_calculator_light(traj);
% % ===== END ensemble and time-averaged MSD calculation =====


save(savefile,'info','traj','msd','-append');

toc
