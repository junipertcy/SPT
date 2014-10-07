d12_b=[];
MSDj_b=[];
j_star_b=[];
MSDj_n_b=[];
MSDjj_b=[];
MSDjjj_b=[];
cont_b=0;
MSDt_b=zeros(1,1);

[ msd_b ] = MSD_calculator_light(traj_b);

for t=1:1
   for j=1:numel(msd_b)
       if and(length(msd_b{j}{1})>99,msd_b{j}{1}(1)>0.00)
           j_star_b = [j_star_b j];
       cont_b=cont_b+1;
       MSDj_b = [MSDj_b sum([msd_b{j}{2}{t} msd_b{j}{3}{t}].^2,2)'];  
       d12_b(cont_b)=(msd_b{j}{1}(2)-msd_b{j}{1}(1))/4/0.01;       
       
%        for tt=1:100
%        MSDjj_b = [MSDjj_b MSDj_b(1:tt)];
%        
%        MSDjjj_b{tt,cont_b} = MSDjj_b;
%        MSDjj_b=[]
%        
%        end
       
       
       else
       
       end
   
   end
    
   
    
   MSDt_b(t) = mean(MSDj_b);
   MSDj_n_b{t} = MSDj_b/MSDt_b(t); 
   
   MSDj=[];
   
end




% MSDjjjj_b=[];
% for j=1:100
% for i=1:cont_b
% MSDjjjj_b = [MSDjjjj_b MSDjjj_b{j,i}];
% 
% end
% Mj(j) = mean(MSDjjjj_b);
% MSDjjjj_b=[];
% end
% 






cont_b=cont_b/t;
j_star_b= j_star_b(1:cont_b);
d12_b = d12_b(1:cont_b);


kk_1_b=[]; for i=j_star_b; kk_1_b = [kk_1_b ; msd_b{i}{2}{1}]; end
kk_10_b=[]; for i=j_star_b; kk_10_b = [kk_10_b ; msd_b{i}{2}{10}]; end
figure(); hold on; [a b] = hist(kk_1_b,100); line(b,a/trapz(b,a));
[a b] = hist(kk_10_b,100); line(b,a/trapz(b,a));


