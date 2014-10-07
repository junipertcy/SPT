% fickian_1.m (Unfinished)
%
% Easy scrpt to access the evolution of the non-gaussian parameter
%
%
% Tzu-Chi Yen, IAMS
% Last modified Sep 2014



% fickian_plot
vvv=[];for i=1:1000;[a b] =hist(displacement{1}{i,1},100); [a1 a2 a3 a4]=histStat(b,a'); vvv(i,:,:,:,:)=[a1 a2 a3 a4];end
vvv2=[];for i=1:1:1000;[a b] =hist(displacement{1}{i,2},100); [a1 a2 a3 a4]=histStat(b,a'); vvv2(i,:,:,:,:)=[a1 a2 a3 a4];end

for i=1:1000;alpaa(i) = vvv(i,4)/3/vvv(i,2)^2-1;end
for i=1:1000;alpaa2(i) = vvv2(i,4)/3/vvv2(i,2)^2-1;end




displacement = displacement_;

displacement_ = displacement;
displacement=[];
% Static frame displacement
for i=1:100
    static_f1{i} = displacement{1}{i,1}(1:round(length(displacement{1}{i,1})*(1/i)));
    static_f2{i} = displacement{1}{i,2}(1:round(length(displacement{1}{i,2})*(1/i)));
    static_f3{i} = displacement{1}{i,3}(1:round(length(displacement{1}{i,3})*(1/i)));
end
displacement = {[static_f2;static_f3]'};