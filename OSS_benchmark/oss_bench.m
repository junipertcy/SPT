% oss_bench.m
%
% Script that plots single-exponential fits to the OSS benchmark experiment.
%
% NOTE:
%      Should only be excuted with traj_* data on the background
%
% Tzu-Chi Yen
% last modified Sep 2014

targ = [];

for w = 1:15 
str_traj = sprintf('%s%d%s','traj_',w,';');
str = ['traj_k = ',str_traj];
eval(str);  
    
str_peaks = sprintf('%s%d%s','peaks_',w,';');
str = ['peaks_k = ',str_peaks];
eval(str);
    
    l_traj = length(traj_k);

    for i=1:l_traj
    ll_traj(i) = length(traj_k{i});
    end

    for j = 1:l_traj

        for i=1:length(peaks_k)
            k=find(peaks_k{i}==traj_k{j}(1,1));
            if isempty(k)==0
                sss(j)=i;
                break
            end
        end
    end


    for i = 1:l_traj
        %targ = [targ  sss(i):sss(i)+ll_traj(i)-1];
        % for density time-evolution
        targ = [targ  sss(i):sss(i)+ll_traj(i)];
    end
end

[a b] = hist(targ,1000);
figure(13);line(b,a/trapz(b,a));hold on;

[k1 k2] = max(a)
b(k2)

a_ = a/trapz(b,a);

% fitting exp
f = fit((b(k2:end)-b(k2))',a_(k2:end)','exp1')
plot(f,'-')
