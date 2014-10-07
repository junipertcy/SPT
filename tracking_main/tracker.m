% tracker.m
%
% Single dye tracking code for fluorescence images, loading
% necessary functions for localization.
%
% User-defined variables are in *run.m* script
%
% Function calls:
%       sif2tif.m
%       detect_particles_IAMS_*.m
%       radialcenter.m
%       saturate.m
%       link_trajectories_IAMS.m   
%       sort_trajectories_IAMS.m
%
%
% revision: Tzu-Chi Yen 11-13-2013
% 1. deleting all unnecessary comments, making code compact 
% 2. adding save functionality for file names with parameters


peaks = [];

%=============================================================
%  Transform *.sif to *.tif
%=============================================================
% write *.tif image stack
for k = N_start:N_start+N_img-1
img(:,:,k) = sif2tif(filename, k, Width, Height, resolution, NumFrames);
%img(:,:,k) = single(imread(filename,'index',k));
%img = img_old(350:500,250:400,:);
%imwrite(uint16(img(:,:,k)),savevideo,'WriteMode','append');
end

for k = N_start:N_start+N_img-1
    
    % NOTE we do not need to normalize the image here
    % Within Gaussian-tracker, an offset is defined within the code
    %maxint(k)=max(max(img(:,:,k)));
    %minint(k)=min(min(img(:,:,k)));
    %images_nor = (img(:,:,k)-minint(k))./(maxint(k)-minint(k));
disp(sprintf('\nParticle recoginition in image %d of %d',k,N_img));
peak = detect_particles_IAMS_radial(img(:,:,k),w,ww,k_rank,pth,pth_int);
peaks = [peaks, peak];
time (k) = k*frame_time;
end

%=============================================================
%                    link trajectories
%=============================================================
viz=0;
peaks = link_trajectories_IAMS(peaks, L, viz, N_end-N_start+1);
[trajectories, all_part] = sort_trajectories_IAMS(peaks);

%=============================================================
% visualize paths  
%=============================================================
figure;
for k=1:N_img-1
    
    oldind = find(peaks{k}(:,6)>0);
    curind = peaks{k}(oldind,6);
    
    X = [peaks{k}(oldind,1), peaks{k+1}(curind,1)].*pixel_size;
    Y = [peaks{k}(oldind,2), peaks{k+1}(curind,2)].*pixel_size;
    
    hand = line(X',Y');
        
    set(hand(:),'Color',[1 0 0]);
    set(hand(:),'LineWidth',[1.0]);
end;

hold off
daspect([1 1 1]);
xlabel('X position (\mum)') 
ylabel('Y position (\mum)') 

set(gca,'fontsize',16)
h_xlabel = get(gca,'xlabel');
set(h_xlabel,'FontSize',20); 
h_ylabel = get(gca,'ylabel');
set(h_ylabel,'FontSize',20); 
box on;


str_peaks = sprintf('%s%d','peaks_',i_spec);
str = [str_peaks, '= peaks;'];
eval(str);

%save(savefile,'peaks','all_part','trajectories','-append');
save(savefile,str_peaks,'-append');

