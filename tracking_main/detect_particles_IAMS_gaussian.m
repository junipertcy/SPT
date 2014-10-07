% detect_particles_IAMS_gaussian.m
% 
% localize sinlge particles using the gaussian-fitting method
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
% 2014/03/28 Tzu-Chi Yen
% Localization kernel modified from detect_particles_IAMS_v2
% 
%
% 2013/07/18 Chia-Lung Hsieh
% Detection of particles above threshold was simplified.
% Use of the function "imhist" is avoided.
% The value for pth needs to be more accurate to avoid multiple detections
%

function peak = detect_particles_IAMS_gaussian(orig,w,ww,k_rank,pth)
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


%  figure(101); hold on;imagesc(orig);colormap(gray);
%  figure(102);hold on; imagesc(filtered);colormap(gray);
%  figure(103); imagesc(dil);colormap(gray);

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
% STEP 3: 2D Gaussian fit
%====================================================================== 

counter = 0;
for ipart=1:npart,
    [n_orig,m_orig]=size(orig);
    if and(and((R(ipart)-ww>0),(C(ipart)-ww>0)),and((R(ipart)+ww<n_orig),(C(ipart)+ww<m_orig)))
    
    counter = counter +1;
    
    % note: LSQCURVEFIT only accepts inputs of data type double
    I_g = double(orig(R(ipart)-ww:R(ipart)+ww,C(ipart)-ww:C(ipart)+ww));
    [n_g,m_g] = size(I_g);
    [X_g,Y_g]=meshgrid(1:n_g,1:m_g);%x-y coordinates
    x_g(:,1)=X_g(:); % x= first column
    x_g(:,2)=Y_g(:); % y= second column
    f_g=I_g(:); % your data f(x,y) (in column vector)

    % Gaussian fitting
    fun = @(c_g,x_g) c_g(1)+c_g(2)*exp(-((x_g(:,1)-c_g(3))/c_g(5)).^2-((x_g(:,2)-c_g(4))/c_g(6)).^2);

    %--- solve with lsqcurvefit
    options=optimset('Display','off','TolX',1e-8,'TolFun',1e-10,'MaxFunEvals',1000);
    
    init_offset = mean([I_g(1:2*ww+1,1)' I_g(1:2*ww+1,2*ww+1)' I_g(1,2:2*ww) I_g(2*ww+1,2:2*ww)]);
    init_peak = max(max(I_g));
    c0_g = [init_offset init_peak ww ww w-1 w-1]; %start-guess here
    
    % MODIFIED: SIGI, 8/11/2012
    [cc_g,resnorm,~,~,~,~,jacobian] = lsqcurvefit(fun,c0_g,x_g,f_g,[],[],options);
    mse = resnorm/(numel(f_g)-length(cc_g));
    covar=jacobian'*jacobian;
    covar=covar\eye(size(covar))*mse;
    % the variable *se* contains the standard deviations for the fitting parameters in the same order
    se=sqrt(diag(covar));     
    % END
    
    Ifit=fun(cc_g,x_g); %your fitted gaussian in vector
    Ifit=reshape(Ifit,[n_g m_g]);%gaussian reshaped as matrix

    %figure; imagesc(I_g);
    %figure; imagesc(Ifit);

    epsx = cc_g(4)-ww-1;
    epsy = cc_g(3)-ww-1;

     R_par(counter) = R(ipart)+epsx;
     C_par(counter) = C(ipart)+epsy;
     I_par(counter) = cc_g(2);      % intensity profile (constrast)
     
     R_se(counter) = se(4); % std for y-position
     C_se(counter) = se(3); % std for x-position
     

    else
    end

end;	


if counter == 0
    
    peak = zeros(1,6);
else
peak = zeros(counter,6);
% checked on 03-28-2014 by TCY
peak(:,1) = C_par;       % row position (x-axis)
peak(:,2) = R_par;       % col position (y-axis)
peak(:,5) = I_par;       % intensity profile

peak(:,3) = C_se;       % std for x-position
peak(:,4) = R_se;       % std for y-position


end
peak = {peak};

return

