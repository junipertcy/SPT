% fickian_2.m (Un-finished)
%
% 1. Script for calculating step size distribution for virus particles on
%    membrane
%
% 2. Non-gaussian parameter can also be accessed in this script
%
% Note:
%       This code loads pre-calculated trajectory data, and then calculates
%       its MSD and displacement distribution.
%
%
%
% Tzu-Chi Yen
% last modified: Jul. 24, 2014
% began: Apr. 22, 2014
%


% ===== Begin User-defined Parameters =====
frame_rate = 100;
n_videos = 38;
str_path = '/Users/junipe/mnt/iams-lustre/tzuchi/SPT_code/input/virus/20140320/vaccinia_dA26_virus/';
str_filename_main = 'vacinnia_dA26_virustrackdata_Cine';

msd = cell(1,n_videos);
displacement = cell(1,n_videos);

tic
for i_spec=1:n_videos
    str_filename = sprintf('%s%d',str_filename_main,i_spec); 
    filename = sprintf('%s%s%s',str_path,str_filename,'.mat');
    load(filename);
    
    
    % Extract single, longest trajectory for analysing
    l_virus = zeros(1,length(trajectories));
    for i=1:length(trajectories)
        l_virus(i) = length(trajectories{i});
    end
    [val idx] = max(l_virus);
    
    traj = trajectories{idx}(1:3000,1:2);
    

    % MSD Calculation
    [msd_all] = MSD_calculator({traj});
    msd{i_spec} = msd_all;

    
    % Displacements calculation
    displacement{i_spec} = traj_displacements(msd{i_spec});
    
    save(str_filename_main,'str_filename_main','displacement','msd');
    
end
toc

% Non-Gaussian parameter
a2 = zeros(1,10);
for j=1:1
    for i=[1:10]
    [ss1, ss2] = hist(displacement{1}{1}.^0.5,10);
    ss1 = ss1./ss2;
    bar([-ss2,ss2],[ss1,ss1])

    a2(j,i) = (moment(ss1,4)/3/moment(ss1,2)^2)-1;
    end
end

%% Plotting
step_slice = [1,10,20];
count_n = 0;
for i = step_slice

    count_n = count_n+1;
    [num_ bin_] = hist(displacement{3}{i},50);
    num_= num_./bin_;
    
    num_norm{count_n} = num_;
    bin_norm{count_n} = bin_;
    
end

figure(2);hold on;
plot(bin_norm{1},num_norm{1},'d','color','k');
plot(bin_norm{2},num_norm{2},'+','color','g');
plot(bin_norm{3},num_norm{3},'o','color','c');
plot(bin_norm{4},num_norm{4},'*','color','b');
plot(bin_norm{5},num_norm{5},'.','color','r');
plot(bin_norm{6},num_norm{6},'o','color','b');
plot(bin_norm{7},num_norm{7},'d','color','k');





x = [0:0.001:0.03];
s = mean_msd{1}(2)-mean_msd{1}(1); 


f = gauss_distribution(x, 0, s);
plot(x,log(f),'.')
grid on
title('Bell Curve')
xlabel('Randomly produced numbers')
ylabel('Gauss Distribution')


ll=[];
for j=1:500
    for i=1:38
    ll= [ll; displacement{i}{j}];
    end
    distp{j}=ll;
    ll=[];
end






