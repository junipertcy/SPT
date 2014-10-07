% vcorr.m
% 
% Script that calculates velocity autocorrelation function to 
% distinguish kinds of subdiffusion 
% (one trajectory outputs a single autocorrelation data)
%
% 
% Refer: S. Burov et al., Phys. Chem. Chem. Phys., 2011, 13, 1800
%        S.C. Weber et al., Biophys. J, 2012, 102, 2443
%
% INPUT: 
%       2*N array of trajectory time-series data (example: traj)
% 
% OUTPUT: 
%       vcor(ti,ep)
%
%
%
% Tzu-Chi Yen, IAMS
% last modified Sep 2014
% began Jan 2014

l_traj = length(traj);

% define delta-time for transient velocity
d_time = 1;

% define transient velocity at certain time
trans_v = zeros(l_traj-d_time,2);
for i=1:l_traj-d_time
    trans_v(i,:) = d_time^-1*(traj(i+d_time,:)-traj(i,:));
end

% (transient) velocity auto-correlation
ll_traj = size(trans_v,1);

for j=0:10
    for k=1:ll_traj-j
    temp_c(k) = dot(trans_v(k+j,:),trans_v(k,:));
    end
    v_corr(j+1) = mean(temp_c);
    temp_c=[];
end

for j=1:100
    for i=1:100
    vcor_tmp(i,j) = dot(traj(j+1+i,:)-traj(j+i,:),traj(i+1,:)-traj(i,:));
    end
end

vcor_1 = l_traj^-1*squeeze(sum(vcor_tmp,1));
vcor_0 = l_traj^-1*squeeze(sum(vcor_tmp(:,1,:),1));
vcor = bsxfun(@rdivide,vcor_1,vcor_0');