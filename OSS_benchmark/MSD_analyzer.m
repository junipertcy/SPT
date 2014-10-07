% MSD_analyzer.m
%
% filter that selects trajectory under certain constraints
% 
% OLD-code, functions may coincide with MSD_collector.m
% used in run.m with minimal constrains in order not to affect the
% final truncation at MSD_collector.m 
% 
%
% Tzu-Chi Yen
% last modified Sep 2014
%
%


% define magic distance
wugaga = 0.0;	% unit: pixel (distance between the maximum
                % trajectory excursion point to that of the mean
                % excursion point.)
min_length=10;  % minimum trajectory length being accepted


% sorting all trajectories from the longest to shortest
[temp traj_sort] = sort(all_part(:,2)-all_part(:,1),'descend');
l_traj_sort = length(traj_sort);

% saving all useful trajectories to one variable
mmm = 0;
for i=1:l_traj_sort
ii = length(trajectories{i})
	if ii>min_length
	mmm = mmm+1;
	
        % First two column: point locations
        traj{mmm} = trajectories{i}(:,1:2);
        % fith column: peak profile
        %traj_I{mmm} = trajectories{i}(:,5);        
    else 
		break
	end
end

lll=[];
% begin wugaga trimming
for i=1:length(traj)
lll(i)=length(traj{i});
end
l_traj_sort = length(lll);

stat_1=[];
for i=1:l_traj_sort
cm_traj{i}=mean(traj{i});       
stat_1(i)=max(sum(bsxfun(@minus,traj{i},cm_traj{i}).^2,2).^0.5);     
end;
                                
cont_n=1;
for i=1:l_traj_sort
    if stat_1(i)>wugaga
        traj_new{cont_n}=traj{i};
        traj_new_I{cont_n}=traj_I{i};
        cont_n=cont_n+1;
    end 
end;

traj=[];
traj_I=[];

traj=traj_new;
traj_new=[];

traj_I=traj_new_I;
traj_new_I=[];

str_traj = sprintf('%s%d','traj_',i_spec);
str = [str_traj, '= traj;'];
eval(str);

str_traj_I = sprintf('%s%d','traj_I_',i_spec);
str = [str_traj_I, '= traj_I;'];
eval(str);






%save(savefile,'sd_all','sd_traj','traj_sort',str_traj,'-append');
save(savefile,str_traj,str_traj_I,'-append');



% short code for truncation
% j=0;for i=1:558
%     if length(trajectories{i}) > 50
%     j=j+1;
%     traj{j}=trajectories{i};
%     end
% end





