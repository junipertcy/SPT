% RWer_loc.m
%
% Script that generates simulated Brownian motion from Gaussian propagators
% The length/time_between_each_jump/estimated_D_constant can be
% adjusted, or loading from a pre-existing variable
%
% NOTE: 
%       output trajectory has the unit of meters
%
% OUTPUT:
%       traj_b: as a cell, similar to variable traj used in other codes
% 
% This script is intended for generating RWer using different propagators 
% (currently only Gaussian is considered)
%
% Tzu-Chi Yen
% last modified Sep. 2014
% starting on 01/20/2014
%
% The backbone of the script is from:
% http://www.advancedlab.org/mediawiki/index.php/Simulating_Brownian_Motion

% % theoretical value of D
% d    = 1.0e-6;              % radius in meters
% eta  = 1.0e-3;              % viscosity of water in SI units (Pascal-seconds) at 293 K
% kB   = 1.38e-23;            % Boltzmann constant
% T    = 293;                 % Temperature in degrees Kelvin
% D    = kB * T / (3 * pi * eta * d)




traj_b=[];

% Using a for loop to generate multiple data sets
for ww=1:cont
    
% ==================================================
% User defined variables
% ==================================================

particleCount = 1;
N = length(msd{j_star(ww)}{1})+1;
dimensions = 2;         % two dimensional simulation
tau = .01;             % time interval in seconds
D = ((msd{j_star(ww)}{1}(2)-msd{j_star(ww)}{1}(1))/4/tau)*1e-12;                % diffusion constant (unit: m^2/sec)
k = sqrt(D*dimensions*tau); % unit: m
time = 0:tau:(N-1)*tau;




    n_count=1;
    particle = { };             % create an empty cell array to hold the results


    for i = 1:particleCount
        particle{i} = struct();
        particle{i}.dx = k*randn(1,N);
        particle{i}.x = cumsum(particle{i}.dx);
        particle{i}.dy = k*randn(1,N);
        particle{i}.y = cumsum(particle{i}.dy);
    end
traj_b{ww} = [particle{i}.x; particle{i}.y]'./(0.096*10^-6);


end


