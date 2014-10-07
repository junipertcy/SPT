% exp_waiting.m
%
% Script that accesses the experimental waiting time distribution by defining
% a one-step distance threshold to the trajectory under study.
%
%
% The waiting time distribution should be used to study CTRW.
% refer: Weigel, A. V. et al, Proc. Natl. Acad. Sci. U.S.A. 108, 6438–6443 (2011).
%
% Inportant OUTPUT:
%                  e_waiting: waiting times at shortest time-lag,
%                  unit: frame time (here: 0.01 sec)
%
%
%
% Tzu-Chi Yen
% last modified Sep 2014
%

% unit: micron
% For a jump size less than the cutoff, it will be considered as
% being trapped; otherwise, it takes a JUMP.
cutoff=0.025;


d12=[];
MSDj=[];
j_star=[];
MSDj_n=[];
e_waiting=[];
cont=0;
MSDt=zeros(1,1);
   for j=1:numel(msd)
       % note: mean(std(traj{j})*0.096) is the mean standard
       % deviation of the trajectory, referring to the OSS
       % experiment, for immobile particles, it should be ~ 0.0345 (micron)
       if and(length(msd{j}{1})>99,mean(std(traj{j})*0.096)>0.0345)
           j_star = [j_star j]; % index that satisfies the
                                % criterion above
       cont=cont+1;
       
       % 
       MSDj = [MSDj sum([msd{j}{2}{1} msd{j}{3}{1}].^2,2)'];  
       
       
       indd = find(MSDj>cutoff);
       MSDj=[];
       
       % waiting times, unit of frame time; in our dye tracking
       % experiment, it is 0.01 sec.
       e_waiting = [e_waiting diff(indd)];
       
       else
       end
   
   end    
   
% here plots the waiting time distribution
figure(1); hold on; [a b]=hist(e_waiting,100);
aa= [log(b);log(a)]';
aa(aa(:,2)<-100,:)=[];
plot(aa(:,1),aa(:,2),'o')