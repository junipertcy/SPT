% detect_particles_IAMS_centroid.m
% 
% localize sinlge particles using the centroid method
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
% 2013/12/30 Tzu-Chi Yen
% Localization kernel modified from detect_particles_IAMS_v2
% The variable Aij modified from I.F. Sbalzarini's original code
%
% 2013/07/18 Chia-Lung Hsieh 
% Detection of particles above threshold was simplified.
% Use of the function "imhist" is avoided.
% The value for pth needs to be more accurate to avoid multiple detections
%

function peak = detect_particles_IAMS_centroid(orig,w,ww,k_rank,pth)

% clear;
% temp = double(imread('Cine06\img_00001.tif'));
% 
% w = 5;         % radius of neighborhood: > particle radius, 
%                % < interparticle distance
% cutoff = 0;  % probability cutoff for non-particle discrimination
% pth = 10;     % upper intensity percentile for detecting particles (rough estimation)
% L = 2;        % maximum displacement between frames
% N_bit = 12; 
% pixel_size = 48e-3; % 48nm, v711, 20um pixel size, 75cn lens, 100x OBJ
% 
% maxint = 3240;
% minint = 3010;
% 
% temp_inv = 2^N_bit - temp;
% orig = (temp_inv-(2^N_bit-maxint))./(maxint-minint);



% some often used quantities
idx = [-w:1:w];     % index vector
dm = 2*w+1;         % diameter
im = repmat(idx',1,dm);
jm = repmat(idx,dm,1);
imjm2 = im.^2+jm.^2;
siz = size(orig);   % image size

%====================================================================== 
% STEP 1: Image restoration
%====================================================================== 

% build kernel K for background extraction and noise removal
% (eq. [4])
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

dil = imdilate(filtered,mask);
filtered=padarray(filtered,[w w]);
dil=padarray(dil,[w w]);


% figure(101); hold on;imagesc(orig);colormap(gray);
% figure(102);hold on; imagesc(filtered);colormap(gray);
% figure(103); imagesc(dil);colormap(gray);

k_rank = k_rank+1;

[Rp,Cp] = find((dil-filtered)==0);
particles = zeros(siz);
max_filtered=max(max(filtered));
suspicious=filtered(sub2ind(siz,Rp,Cp))./max_filtered;
sorted_suspicious = sort(suspicious,'descend');
k_max = sorted_suspicious(min(k_rank,end));
V1 = find(suspicious > k_max);
V2 = find(suspicious > pth);
V=intersect(V1,V2);
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

% zeroth order intensity moment of all particles
m0 = zeros(npart,1);

counter = 0;

% for each particle, compute zeroth order moment and position corrections
% epsx, epsy
for ipart=1:npart
    counter = counter +1;
    width_g = ww; % crop a small image around the particle for fitting 
    epsx = 1; epsy = 1;
    while or(abs(epsx)>0.5,abs(epsy)>0.5),
	% lower and upper index bounds for all particle neighborhoods
	% in local coordinates. Recalculate after every change in R,C
	li = 1-(R-width_g-saturate(R-width_g,1,siz(1)));
	lj = 1-(C-width_g-saturate(C-width_g,1,siz(2)));
	ui = dm-(R+width_g-saturate(R+width_g,1,siz(1)));
	uj = dm-(C+width_g-saturate(C+width_g,1,siz(2)));
	% masked image part containing the particle 
        % (original code uses filtered, which seems wrong)
    Aij = orig(R(ipart)+li(ipart)-width_g-1:R(ipart)+ui(ipart)-width_g-1,...
	    C(ipart)+lj(ipart)-width_g-1:C(ipart)+uj(ipart)-width_g-1).* ...
	    mask(li(ipart):ui(ipart),lj(ipart):uj(ipart));
	% moments
	m0(ipart) = sum(sum(Aij));    % eq. [6]
	% eq. [7]
	m2(ipart) = sum(sum(imjm2(li(ipart):ui(ipart),lj(ipart):uj(ipart))...
	    .*Aij))/m0(ipart); 
	% position correction
	epsx = sum(sum(im(li(ipart):ui(ipart),lj(ipart):uj(ipart)).*Aij))/m0(ipart);
	epsy = sum(idx(lj(ipart):uj(ipart)).*sum(Aij))/m0(ipart);
	% if correction is > 0.5, move candidate location
	if abs(epsx)>0.5,
	    R(ipart) = R(ipart)+sign(epsx);
	end;
	if abs(epsy)>0.5,
	    C(ipart) = C(ipart)+sign(epsy);
	end;
    end; 
    
    R_par(counter) = R(ipart)+epsx;
    C_par(counter) = C(ipart)+epsy;
    
    % Intensity profile
    % peak intensity - background (assuming ww is large enough)
    siz_A = size(Aij);
    www_c = siz_A(1);
    www_r = siz_A(2);
    cycle_A = [Aij(1:www_c,1)' Aij(1:www_c,www_r)' Aij(1,2:www_r-1) Aij(www_c,2:www_r-1)];
    I_par(counter) = max(max(Aij))-mean(cycle_A);

    
end 


if counter == 0   
    peak = zeros(1,2);
else
peak = zeros(counter,2);
peak(:,1) = C_par;       % row position
peak(:,2) = R_par;       % col position
peak(:,5) = I_par;       % intensity

end


peak = {peak};


return

