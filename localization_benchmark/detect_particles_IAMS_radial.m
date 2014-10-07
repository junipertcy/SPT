% detect_particles_IAMS_radial.m
% 
% localize sinlge particles using the radial symmetry method (OLD)
%
% INPUT:
%       orig  : input image
%       w     : size for gaussian convolution (default: 3)
%       ww    : cropped size for localization (default: 6)
%       k_rank: threshold for number of particles detected        [>1]
%       pth   : threshold for intensity of the particles detected [0-1]
%
% OUTPUT:
%       peak: cell array for detected particles of in x- and y- coordinates
%
%
% 03-27-2014 Tzu-Chi Yen
% Last modification. Checking functionality.
%
% 01-13-2014 Tzu-Chi Yen
% Kernel change: implementing detection scheme based on radial symmetry,
% developed by R.P.
%
% 2013/07/18 Chia-Lung Hsieh 
% Detection of particles above threshold was simplified.
% Use of the function "imhist" is avoided.
% The value for pth needs to be more accurate to avoid multiple detections

function [peak] = detect_particles_IAMS_radial(orig,w,ww,k_rank,pth)

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
% refer: I.F. Sbalzarini et al, J. of Structural Biology 151 (2005) 182-195
% build kernel K for background extraction and noise removal
% (eq. [4])
lambdan = 1;    % digitalization noise in the CCD camara and the frame grabber (why close to 1 pixel??)
B = sum(exp(-(idx.^2/(4*lambdan^2))));
B = B^2;
K0 = 1/B*sum(exp(-(idx.^2/(2*lambdan^2))))^2-(B/(dm^2));
K = (exp(-(imjm2/(4*lambdan^2)))/B-(1/(dm^2)))/K0;
figure(31); imagesc(K);
% apply convolution filter
filtered = conv2(orig,K,'valid');


%====================================================================== 
% STEP 2: Locating particles
%====================================================================== 

%[R C] = find (filtered(:,:)> 1-pth*0.01);

% % determining upper pth-th percentile of intensity values
% pth = 0.01*pth;
% 
% sort_filtered = sort(filtered(:));
% l2 = length(sort_filtered);
% thresh = sort_filtered(round(l2*pth));

% [cnts,bins] = imhist(filtered);
% l = length(cnts);
% k = 1;
% while sum(cnts(l-k:l))/sum(cnts) < pth,
%     k = k + 1;
% end;
% thresh = bins(l-k+1);

% generate circular mask of radius w
mask = zeros(dm,dm);
mask(find(imjm2 <= w*w)) = 1;

% identify individual particles as local maxima in a
% w-neighborhood that are larger than thresh
dil = imdilate(filtered,mask);
filtered=padarray(filtered,[w w]);
dil=padarray(dil,[w w]);

figure(32); imagesc(filtered);
figure(33); imagesc(dil);


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
%    I_g = orig(R(ipart)-width_g:R(ipart)+width_g,C(ipart)-width_g:C(ipart)+width_g);

% Note the cropped image is different from what was normally used. 
% This is because in the radialcenter.m code, the y-axis inceases downward.
    % [Working within image coordinate, 03-27-2014]
    I_g = orig(C(ipart)-width_g:C(ipart)+width_g,R(ipart)-width_g:R(ipart)+width_g);



    [R_star C_star]=radialcenter(I_g);
    R_par(counter) = R(ipart)-width_g+R_star-1;
    C_par(counter) = C(ipart)-width_g+C_star-1;
    
%     R_par(counter) = R(ipart);
% C_par(counter) = C(ipart);
%     R_par(counter) = R_star;
%     C_par(counter) = C_star;

    
    
    else
        
    end

end;



if counter == 0   
    peak = zeros(1,2);
else
peak = zeros(counter,2);
peak(:,1) = R_par;       % row position
peak(:,2) = C_par;       % col position


end
peak = {peak};
return

