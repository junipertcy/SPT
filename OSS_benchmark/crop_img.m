% crop_img.m
%
% This MATLAB script reads location data in *ttt*, and crop
% partial video clips centered at each of the specified location points.
%
%
% NOTE: 
%      im_loc needs to be excuted for blinking sites to be detected
%
%
%
% 
% Tzu-Chi Yen
% last modified Sep 2014
% began Jul 2014

tic


ikm_total=[];


% ===== Begin User-defined Parameters =====
for w = 1   % defining the name of variable to be analysed; say,
             % if w=13, then the N*2 trajectory variable traj_13
             % must be pre-exist in the memory
             % w can be a vector, say, 1:15
    
str_traj = sprintf('%s%d%s','traj_',w,';');
str = ['traj_k = ',str_traj];
eval(str); 


% the mean location of each trajectory
mmm=[];
for i=1:length(traj_k)
mmm(:,:,i)=mean(traj_k{i});
end
mmm = squeeze(mmm)';

[tart] = im_loc(mmm,1);

%str_path = '/Users/junipe/mnt/test_images/';
str_path = '/Users/junipe/Desktop/';
str_filename_main = 'PCA_PCD_Trolox_';
n_videos = w;      % number of videos analysed

% cropped image clip size = (2*w+1)-square
w_ = 1;       	 % radius of neighborhood > particle radius

%ttt=mmm;
ttt=tart;

% Start frame for localization
N_start = 1;
% End frame for localization
N_end = 2000;




% ===== End User-defined Parameters =====

%=============================================================
% Begin SPT calculation
%=============================================================
% First and foremost, the videos will be loaded and fed into a
% localization algorithm to extract the respective trajectories

I_g=[];
for i_spec=n_videos
    str_filename = sprintf('%s%d',str_filename_main,i_spec);    
    filename = sprintf('%s%s%s',str_path,str_filename,'.sif');
info = func_ext_uli_sifread2(filename);
Width = info.Width; 
Height = info.Height; 
resolution = info.resolution;
NumFrames = info.NumFrames;  



N_img = N_end-N_start+1;


    for k = N_start:N_start+N_img-1
    img(:,:,k) = sif2tif(filename, k, Width, Height, resolution, NumFrames);
    end

    for k = N_start:N_start+N_img-1
        ttt=round(ttt);
        
        
        
        
        kk = size(ttt,1);
        orig = img(:,:,k);
        [n_orig,m_orig]=size(orig);
        cont = 0;
        for j = 1:kk
            if and(and((ttt(j,1)-w_>0),(ttt(j,2)-w_>0)),and((ttt(j,1)+w_<n_orig),(ttt(j,2)+ww<m_orig)))
            cont = cont+1;
            I_g(:,:,k,cont)=orig(ttt(j,2)-w_:ttt(j,2)+w_,ttt(j,1)-w_:ttt(j,1)+w_);
            end
        end
    end
end


ikm=[];
for j=1:size(ttt,1);
    for i=1:N_end
        ikm(i,j) = mean2(I_g(:,:,i,j));
        %ikm(i,j) = max(max(I_g(:,:,i,j)));
    end
end

%img=[];
%I_g=[];
traj_k=[];

ikm_total = [ikm_total ikm];

end


figure();plot(ikm,'-')
%=============================================================
% Code for saving cropped images into a tif-stack
%=============================================================

% for i=1:N_end
% imwrite(uint16(I_g(:,:,i,1)),'PCA_PCD_Trolox_ex_long_wo_blinking_big_1.tif','WriteMode','append');
% %imwrite(uint16(img(:,:,i)),'savevideo3.tif','WriteMode','append');
% 
% end
% % 
toc

% traj=traj_;
% 
% m=0;
% for i=1:24299
% if length(msd{i}{1})>99;
%     m=m+1;
%     msd_l{m} = msd{i};
% end
% end
% 
% figure();
% for i=1:100
%    
%     plot(log([1:100]*0.01),log(msd_l{i}{1}(1:100)),'-');
%    hold on;
% end
    