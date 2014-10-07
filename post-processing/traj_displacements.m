% traj_displacements.m
% This function calculates the (unnormalized) displacement distribution 
% at each time-lag for a trajectory set.
%
% INPUT: 
%      msd_all: cell array from MSD_calculator.m
%
% OUTPUT:
%      displacement:
%
% Tzu-Chi Yen
% Last modified: 04-28-2014
% Feb. 17, 2014


function [ displacement ] = traj_displacements( msd_all, tlag)

% time-lag for calculating displacements
%tlag=10000;
%tlag = 20;
    
% pre-allocate memory    
displacement = cell(tlag,3);
    
% displacement distribution
disp_l=[];disp_x=[];disp_y=[];
    for ii=1:tlag
        
        for i=1:length(msd_all)
            if length(msd_all{i}{3}) >= tlag
                disp_l = [disp_l; msd_all{i}{3}{ii}];
                disp_x = [disp_x; msd_all{i}{4}{ii}];
                disp_y = [disp_y; msd_all{i}{5}{ii}];
            else
                error(['Error: time window is setting too long, ' ...
                       'trajectories are not long enough!'])
            end
        end
        displacement{ii,1}=disp_l;
        displacement{ii,2}=disp_x;
        displacement{ii,3}=disp_y;
        disp_l=[];disp_x=[];disp_y=[];
    end

end
