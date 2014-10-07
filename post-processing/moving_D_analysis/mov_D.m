% mov_D.m
%
% Script for 1. Dividing trajectories into sub-trajectories
%            2. Calculate the diffusivity therein
%            3. Comparison with Brownian motion
%
% INPUT:
%       traj: trajectory
%       traj_b: Brownian trajectory 
%
% OUTPUT:
%       d12: transient D12 
%       err: D12-error
%       d12_b: transient D12 for the brownian motion 
%       err_b: D12-error for the brownian motion
%
%
%
%   Tzu-Chi Yen
%   Last modified Sep 2014
%   Project began at Jul. 29, 2014


% load the trajectory
%load('/lustre/lwork/yhlin/tzuchi/SPT_code/output/vacinnia/20140717/r_sum_dopc_5_[4_4_3_0.3_10_10000].mat');
%msd = MSD_calculator_light({traj_10000{1}(1:2000,:)});
%RWer_loc;

traj1 = traj(:,1);
traj2 = traj(:,2);
traj1_b = traj_b(1:1000,1);
traj2_b = traj_b(1:1000,2);




%% user defined parameters
cs = 10;   % chunksize
sh = 1;    % shifting size
           % Note: Ying-Hsiu may have experience in this.

%% Initialization
a1=[];
a2=[];

a1_b=[];
a2_b=[];

d12=[];
err=[];

d12_b=[];
err_b=[];

%% Viral Trajectory

a1 = traj1(bsxfun(@plus,(1:cs),(0:sh:length(traj1)-cs)'));
a2 = traj2(bsxfun(@plus,(1:cs),(0:sh:length(traj2)-cs)'));


%% Brownian

a1_b = traj1_b(bsxfun(@plus,(1:cs),(0:sh:length(traj1_b)-cs)'));
a2_b = traj2_b(bsxfun(@plus,(1:cs),(0:sh:length(traj2_b)-cs)'));


for i = 1 : size(a1,1)
traj_t = [a1(i,:); a2(i,:)]';
msd = MSD_calculator_light({traj_t});
d12(i) = (msd{1}{1}(2)-msd{1}{1}(1))/4/0.0005;
err(i) = d12(i)*(3/(cs-2))^0.5;
i
traj_tb = [a1_b(i,:); a2_b(i,:)]';
msd = MSD_calculator_light({traj_tb});
d12_b(i) = (msd{1}{1}(2)-msd{1}{1}(1))/4/0.0005;
err_b(i) = d12_b(i)*(3/(cs-2))^0.5;
end

%=============================================================
% parameters for saving the output
%=============================================================
%str_cs = num2str(cs);
%savefile = ['/lustre/lwork/yhlin/tzuchi/SPT_code/output/vacinnia/20140717/','mov_D_new_',str_cs,'.mat'];

%save(savefile,'d12','err','d12_b','err_b','traj_t','traj_tb');



%scatterError([1:size(a1,1)]*cs,d12,zeros(size(a1,1),1),err); hold on;