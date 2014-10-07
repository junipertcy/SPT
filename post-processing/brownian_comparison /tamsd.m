% tamsd.m
%
% 1. Script that calculates the time-averaged mean square displacement
% of a single trajectory
% 2. Calculates the D-alpha scatter plot (w/ or w/o errorbars)
%
%
% NOTE:
%      should be excuted after MSD_collector.m has finished
%      trajectory and msd processing.
%
% Tzu-Chi Yen, IAMS
% last modified Sep 2014
%

m=0;
p=[];
%i=find(and(traj_var>5, traj_var<100))
figure();
for i=j_star
    m=m+1;
    hold on; 
    plot(log([1:100]*0.01),log(msd{i}{1}(1:100)),'-');
    x = log([1:20]*0.01);
    y = log(msd{i}{1}(1:20));
    [p(m,:),ErrorEst] = polyfit(x,y,1);
    D_sigma(m) = (3/(length(traj{i})-2))^0.5;
end
plot(log([1:100]*0.01),log(MSDt),'o-');
plot(log([1:20]*0.01),log(MSDt(1:20)),'o-');

%% D-alpha scatter plot
% check with the reduced localization in Qian or Michalet's paper
% (example: Phys Rev E, 2010 vol. 82 (4) p. 041914) 
% to see how many points in the MSD should be used to calculate the
% diffusion constant
D_sigma = D_sigma(1:cont);
scatterError(d12,p(:,1),D_sigma',zeros(numel(j_star),1))
hold on; plot(d12,p(:,1),'o')
