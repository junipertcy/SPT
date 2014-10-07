% conditional.m (Unclear interpretation, and unfinished project)
%
% Scripts that show the correlation of jump size for different time-steps
% Comparison is made with Brownian motion
%
% refer: Guan, J., et al, ACS Nano (2014). doi:10.1021/nn405476t
% for more ideas related to this
%
%
% Tzu-Chi Yen, IAMS
% last modified Sep 2014


RWer_loc;

cs = 10;   % chunksize
sh = 10;

for jj=1:1
for k=1:10
    d1=[];d2=[];d1_b=[];d2_b=[];

    
    traj_t = traj{jj};
    
    
if k<=5
figure(1); subplot(5,1,k);
elseif and(k>5,k<=10)
figure(2); subplot(5,1,k-5);
else
figure(3); subplot(5,1,k-10);
end


for i=1:length(traj_t)-k
    d1(i) = sqrt(sum((traj_t(i,:)-traj_t(i+1,:)).^2));
    d2(i) = sqrt(sum((traj_t(i,:)-traj_t(i+k,:)).^2));

   d1_b(i) = sqrt(sum((traj_b(i,:)-traj_b(i+1,:)).^2));
   d2_b(i) = sqrt(sum((traj_b(i,:)-traj_b(i+k,:)).^2));

end





dd = [d1; d2];
[aa bb]=sort(d1);
dd1=dd(:,bb);
dd_p1 = dd1(1,:);
dd_p2 = dd1(2,:);

corr_d(jj,k) = corr(d1',d2');

a1 = dd_p1(bsxfun(@plus,(1:cs),(0:sh:length(dd_p1)-cs)'));
a2 = dd_p2(bsxfun(@plus,(1:cs),(0:sh:length(dd_p2)-cs)'));
d_res = [mean(a1'); mean(a2')]';


dd_b = [d1_b; d2_b];
[aa_b bb_b]=sort(d1_b);
dd1_b=dd_b(:,bb_b);
dd_p1_b = dd1_b(1,:);
dd_p2_b = dd1_b(2,:);

corr_d_b(k) = corr(d1_b',d2_b');

a1_b = dd_p1_b(bsxfun(@plus,(1:cs),(0:sh:length(dd_p1_b)-cs)'));
a2_b = dd_p2_b(bsxfun(@plus,(1:cs),(0:sh:length(dd_p2_b)-cs)'));
d_res_b = [mean(a1_b'); mean(a2_b')]';




plot(d_res(:,1),d_res(:,2),'o-');
hold on;
plot(d_res_b(:,1),d_res_b(:,2),'x-','color','r');





% figure();scatter(d1,d2)
% hold on;
% scatter(d1_b,d2_b)

end

end


