% acs_SNr.m
%
% Code to access the signal-to-noise ratio of an image
%
%
% INPUT: 
%	  orig: an image as an matlab array
%         w: size for gaussian convolution (default: 2)
%         ww: cropped size for localization (default: 4)
%         k_rank: threshold for number of particles detected     [>1]
%         pth: threshold for intensity of the particles detected [0-1]
%
% OUT:
%		snr: SNr
%
%	
% Tzu-Chi Yen
% last modified 02-11-2014
%
%

function [snr] = acs_SNr(orig,w,ww,k_rank,pth)

% The main code is cropped from usual detect_particles_IAMS_*.m routine

% some often used quantities
idx = [-w:1:w];     % index vector
dm = 2*w+1;         % diameter (for convoluted, filtered image -> Gaussian mask)
im = repmat(idx',1,dm);
jm = repmat(idx,dm,1);
imjm2 = im.^2+jm.^2;
siz = size(orig);   % image size

%====================================================================== 
% STEP 1: Image restoration
%====================================================================== 
lambdan = 1;    % digitalization noise in the CCD camara and the frame grabber (why close to 1 pixel??)
B = sum(exp(-(idx.^2/(4*lambdan^2))));
B = B^2;
K0 = 1/B*sum(exp(-(idx.^2/(2*lambdan^2))))^2-(B/(dm^2));
K = (exp(-(imjm2/(4*lambdan^2)))/B-(1/(dm^2)))/K0;
%figure(31); imagesc(K);
% apply convolution filter
filtered = conv2(orig,K,'valid');

%====================================================================== 
% STEP 2: Locating particles
%====================================================================== 
% generate circular mask of radius w
mask = zeros(dm,dm);
mask(find(imjm2 <= w*w)) = 1;

% identify individual particles as local maxima in a
% w-neighborhood that are larger than thresh
dil = imdilate(filtered,mask);
filtered=padarray(filtered,[w w]);
dil=padarray(dil,[w w]);

%figure(32); imagesc(filtered);
%figure(33); imagesc(dil);


[Rp,Cp] = find((dil-filtered)==0);

particles = zeros(siz);
max_filtered=max(max(filtered));

suspicious=filtered(sub2ind(siz,Rp,Cp))./max_filtered;
sorted_suspicious = sort(suspicious,'descend');
k_max = sorted_suspicious(min(k_rank,end));

V1 = find(suspicious > k_max);
V2 = find(suspicious > pth);
V=intersect(V1,V2); % linear index for the vector *suspicious*
R = Rp(V);
C = Cp(V);
    
%======================================================================
% STEP 2b: Delete "double detection"
%======================================================================

iii = 1;
while iii == 1
     iii = 0;
     kkk = 1;
     for jjj = 2:length(R)
         kkk = kkk + 1;
         if ( abs(R(kkk-1)-R(kkk))+abs(C(kkk-1)-C(kkk)) ) <= w % if two det. pixels are next neighbors (including diagonal)
             iii = 1;
             R(kkk) = [];
             C(kkk) = [];
             kkk = kkk - 1; % to account for now shorter R and C vectors
         end
     end
end

particles(sub2ind(siz,R,C)) = 1;
npart     = length(R);


%====================================================================== 
% STEP 3: Refining location estimates
%====================================================================== 
counter = 0;
for ipart=1:npart
    width_g = ww; % crop a small image around the particle for fitting 
    [n_orig,m_orig]=size(orig);
    if and(and((R(ipart)-width_g>0),(C(ipart)-width_g>0)),and((R(ipart)+width_g<n_orig),(C(ipart)+width_g<m_orig)))
    counter = counter +1;
    I_g = orig(R(ipart)-width_g:R(ipart)+width_g,C(ipart)-width_g:C(ipart)+width_g);
    orig_tmp = orig;
    orig_tmp(R(ipart)-width_g:R(ipart)+width_g,C(ipart)-width_g:C(ipart)+width_g)=0;
    cropped(ipart,:,:) = I_g;  
    else    
    end
end


if counter == 0   
    disp('No peak detected!');
    snr = 0;
else
    
    % calcualte background RMS
    orig_tmp = reshape(orig_tmp,1,siz(1)^2);
    orig_tmp(orig_tmp==0)=[];
%    im_RMS=sqrt(mean(orig_tmp.^2));
    im_RMS=mean(orig_tmp);

    % calculate signal intensity
    if exist('cropped')~=0
%    im_I=mean(max(max(cropped)));
    im_I=mean(max(max(cropped)))-im_RMS;
    else 
        im_I = 0;
    end
%snr = im_I/im_RMS;
snr = im_I^0.5;

end


return



