% MSD_collector.m
%
% filter that selects trajectories under certain constraints
% should be used after MSD_calculator_light is excuted and
% msd-cell-data is obtained
%
%
%
% Tzu-Chi Yen
% Last modified Sep 2014
% Began Aug 2014
%
%


d12=[];MSDj=[];j_star=[];
MSDj_n=[];MSDjj=[];MSDjjj=[];
cont=0;
MSDt=zeros(1,1);
for t=1:50
    
   for j=1:numel(msd)
       % note: mean(std(traj{j})*0.096) is the mean standard
       % deviation of the trajectory, referring to the OSS
       % experiment, for immobile particles, it should be ~ 0.0345 (micron)
       if and(length(msd{j}{1})>49, mean(std(traj{j})*0.096)>0.0345)
           j_star = [j_star j]; % index that satisfies the
                                % criterion above
       cont=cont+1;
       MSDj = [MSDj sum([msd{j}{2}{t} msd{j}{3}{t}].^2,2)'];  
       d12(cont)=(msd{j}{1}(2)-msd{j}{1}(1))/4/0.01;       
       
%        for tt=1:100
%        MSDjj = [MSDjj MSDj(1:tt)];
%        
%        MSDjjj{tt,cont} = MSDjj;
%        MSDjj=[];
%        
%        end
       
       
       else
       
       end
   
   end
    
   
   
   MSDt(t) = mean(MSDj);
   %MSDj_n{t} = MSDj/MSDt(t); 
   
   MSDj=[];
   
end




% MSDjjjj=[];
% for j=1:100
% for i=1:cont
% MSDjjjj = [MSDjjjj MSDjjj{j,i}];
% 
% end
% Mj(j) = mean(MSDjjjj);
% MSDjjjj=[];
% end
for i=1:numel(msd)
    ll(i)=length(traj{i});
end

cont=cont/t;
j_star= j_star(1:cont);
d12 = d12(1:cont);

% Here, displacement distribution is calculated
kk_1=[]; for i=j_star; kk_1 = [kk_1 ; msd{i}{2}{1}]; end
kk_10=[]; for i=j_star; kk_10 = [kk_10 ; msd{i}{2}{10}]; end
figure(); hold on; [a b] = hist(kk_1,100); line(b,a/trapz(b,a));
[a b] = hist(kk_10,100); line(b,a/trapz(b,a));


