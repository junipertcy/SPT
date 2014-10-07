% MSD_calculator_light.m
%
% LIGHT-version: for long trajectory use
% This script calculates the MSD (at certain time-lag),
% displacements x and y, using the stationary frame definition
%
% NOTE: 
%      Specify correct pixel size before usage
%
% INPUT:
%      traj: trajectory data, N*2-vector
%
% OUTPUT:
%      msd_all: an 1*3 cell, MSD, x-displacements, y-displacements
%      at different time-lag; unit: micron^2, micron, micron
%
%
% Tzu-Chi Yen   
%    
% Last modified: 09-10-2014
% This script only serves to accept single trajectory video.



function [ msd_all ] = MSD_calculator_light(traj)
pixel_size = 0.096; % unit: micron, FL tracking experiment
%pixel_size = 0.048; % unit: micron, iSCAT experiment

n_data=length(traj);          % number of trajectories analyzed
msd_all=cell(1,n_data);


for nmm=1:n_data
    XY = traj{nmm}*pixel_size;
    N_img = length(XY(:,1));            % loading trajectory x-y data
    msd_res=cell(1,N_img-1);
    disp_x = cell(1,N_img-1);
    disp_y = cell(1,N_img-1);
    msd_mean=zeros(1,N_img-1);
    MSD_result_x=cell(1,1/2*(N_img-1)*N_img);
    MSD_result_y=cell(1,1/2*(N_img-1)*N_img);
    MSD_result=cell(1,1/2*(N_img-1)*N_img);
    
    for t_step = 1:N_img-1                          % t_step for varying time step (for each trajectory: larger dN, fewer statistical samples)
            index = 1/2*(t_step-1)*t_step + 1;      % lower triangular index number
            MSD_result_x{1,index} = diff(XY(1:t_step:N_img,1));
            MSD_result_y{1,index} = diff(XY(1:t_step:N_img,2));
            MSD_result{1,index} = MSD_result_x{1,index}.^2 + MSD_result_y{1,index}.^2;
    end

 
    % Computation Expensive Part: Need to re-design...
    for i=1:N_img-1
        %tmp = zeros(N_img-i,1);
        tmp=[];
        tmp_x=[];
        tmp_y=[];
        
        for j=1:i
        tmp = [tmp; MSD_result{1/2*(i-1)*i+j}];
        tmp_x = [tmp_x; MSD_result_x{1/2*(i-1)*i+j}];
        tmp_y = [tmp_y; MSD_result_y{1/2*(i-1)*i+j}];
        end
        msd_mean(i)=mean(tmp);
        disp_x{i} = tmp_x;
        disp_y{i} = tmp_y;
    end
    
    % Free the memory
    tmp=[];
    tmp_x=[];
    tmp_y=[];
    
    msd_all{1,nmm}={msd_mean,disp_x,disp_y};      % for analyzing MSD progression over time-lag
    %NOTE: msd_res can be used to plot the CDF and be looked for different diffusion modes                                              
end

end

