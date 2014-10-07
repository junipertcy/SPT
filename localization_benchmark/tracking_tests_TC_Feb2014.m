% tracking_tests_TC_Feb2014.m
%
% Characterization of particle fitting algorithms: In contrast to the
% routine by RP's code, this function is written to benchmark the localization
% machinary used in Chia-Lung Hsieh's lab.
%
% Comparing:
% -- (1) radial symmetry fitting: minimizes distance to gradient line method
% -- (2) centroid calculation
% -- (3) Gaussian fitting via nonlinear least squares
%
% uses
%  fitline.m (RP), psf2d.m, modelimage.m [revised, fast version], in addition to
%  radialcenter.m and other main tracking functions (detect_particles_IAMS_radial.m, 
%           detect_particles_IAMS_centroid.m, detect_particles_IAMS_gauss_nonlin.m)
%
% Inputs
%    Ntrials : number of tests per signal/noise (default 1000)
%    maxx0 : max. shift of true particle image center in each dimention, px (default 0.5)
%    SNrparams : Signal-to-noise ratios to test; either 
%        -- a single number (default; value = 20)
%        -- a 3-element array containing min(SNr), max(SNr), and the number of
%           SNr to examine.
%        -- an array that is not 3-elements in size that contains the 
%           SNr values to examine
%    outfilename : name of MAT file to output (default
%        'track_test_output.mat')
%    w, ww, k_rank, pth: common input parameters for function detect_particles_IAMS_*.m
%          
%
% Outputs
%    sigma: Precision (px) (error minus linear fit)
%    time : execution time per fit (s)
%    bias : bias (slope of error vs. true position)
%    sigbias  : uncertainty in bias
%    toterror : total error of fit (px)
%
% Tzu-Chi Yen
% based on tracking_tests_28Aug2011.m (original version)
%      and tracking_tests_RP_Apr2012.m (RP's modification)
%
% Feb. 6, 2014
% last modified Feb. 6, 2014

function [sigma time bias sigbias toterror] = tracking_tests_TC_Feb2014(Ntrials, maxx0, SNrparams, w,ww, k_rank, pth, outfilename)

% Initialize random number generator
s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);

% Parameters for the simulated images
scale = 96; % CCD pixel scale (nm/pixel)
lambda = 553; % wavelength of light (nm) [the emission peak of atto532 is 553 nm]
NA = 1.3;  % numerical aperture
N = 32; % size of the simulated image (NxN pixels).
bk = 400.0;  % background (dark) intensity.
dhr = 1.0;  % high resolution grid size for simulated images, nm
fs = sprintf('Using dhr = %.1f', dhr); disp(fs)

% Parameters for the tests
if ~exist('Ntrials', 'var') || isempty(Ntrials)
    Ntrials = 1000;  % number of tests
    fs = sprintf('Ntrials = %d',Ntrials); disp(fs)
end
if ~exist('maxx0', 'var') || isempty(maxx0)
    maxx0 = 0.5;  % max. shift of true particle image center, px
    fs = sprintf('maxx0 = %.1f',maxx0); disp(fs)
end
if ~exist('SNrparams', 'var') || isempty(SNrparams)
    SNrparams = 20;  % SNr to examine
    fs = sprintf('Single SNr: %.2f',SNrparams); disp(fs)
end
if ~exist('outfilename', 'var') || isempty(outfilename)
    outfilename = 'track_test_output.mat';  % output file name
end

% Pre-calculate the point spread function (over a "large" area)
maxsize = ((N-1) + 2.0*maxx0)*scale;  % maximum needed size for PSF, nm
npsf = round(maxsize/dhr/2.0);
bigpsf = psf2d(npsf, dhr/1000, NA, 1.3, lambda/1000);  % the point-spread function

if length(SNrparams)==1
    SNrarray = SNrparams;  % Just one value
elseif length(SNrparams)==3
    % the range is specified
    SNrmin = SNrparams(1);  
    SNrmax = SNrparams(2);  
    Ntests = SNrparams(3);  % 20 % Number of signal-noise ratios to examine
    SNrarray = logspace(log10(SNrmin), log10(SNrmax), Ntests);  % All signal-to-noise ratios to examine
else
    % the array is specified
    SNrarray = SNrparams;
end

Nmethods = 3;  % number of methods to test

% Allocate memory
sigma  = zeros(Nmethods, length(SNrarray));
time = zeros(Nmethods, length(SNrarray));
bias = zeros(Nmethods, length(SNrarray));
sigbias = zeros(Nmethods, length(SNrarray));
toterror = zeros(Nmethods, length(SNrarray));

for j=1:length(SNrarray)
    fs = sprintf('j= %d of %d; Signal/noise %.3f', j, length(SNrarray), SNrarray(j)); disp(fs)
    
    disp('Creating model images...');
    % ... each with random centers, Poisson noise
    x0 = 2*maxx0*rand(Ntrials,1)-maxx0;
    y0 = 2*maxx0*rand(Ntrials,1)-maxx0;
    %x0 = maxx0*rand(Ntrials,1);  % one quadrant only
    %y0 = maxx0*rand(Ntrials,1);  % one quadrant only
    %y0 = x0;  % for testing "displacements" at 45 degree angle to x,y axes
    z = modelimage(SNrarray(j), N, x0*scale, y0*scale, bk, scale, lambda, NA, dhr, maxx0*scale, bigpsf);
    
    disp('1 Starting radial symmetry fitting')
    x0_r = zeros(1,Ntrials); y0_r = zeros(1,Ntrials);
    tic
    for k=1:Ntrials
        peak = detect_particles_IAMS_radial(z(:,:,k),w,ww,k_rank,pth);
        x0_r(k) = peak{1}(1);
        y0_r(k) = peak{1}(2);
    end
    time(1,j) = toc;

        
    disp('2 Starting centroid fitting')
    tic
    x0_c = zeros(1,Ntrials); y0_c = zeros(1,Ntrials);
    for k=1:Ntrials
        peak = detect_particles_IAMS_centroid(z(:,:,k),w,ww,k_rank,pth);
        x0_c(k) = peak{1}(1);
        y0_c(k) = peak{1}(2);
    end
    time(2,j) = toc;

%     disp('3 Starting Gaussian fitting via least-squares')
%     x0_g = zeros(1,Ntrials);     y0_g = zeros(1,Ntrials);
%     tic
%     for k=1:Ntrials        [
%         peak = detect_particles_IAMS_detect_particles_IAMS_gauss_nonlin(z(:,:,k),w,k_rank,pth);
%         [x0_g(k) y0_g(k)] = [peak{1}(1) peak{1}(2)];
%     end
%     time(3,j) = toc;

    
    % Characterize the various methods (needlessly redundant code...)
    [bias(1,j) sigbias(1,j) sigma(1,j) toterror(1,j)] = calcbias(x0, y0, x0_r, y0_r, N);
    [bias(2,j) sigbias(2,j) sigma(2,j) toterror(2,j)] = calcbias(x0, y0, x0_c, y0_c, N);
%     [bias(3,j) sigbias(3,j) sigma(3,j) toterror(3,j)] = calcbias(x0, y0, x0_g, y0_g, N);

    % Time per fit
    time(:,j) = time(:,j)/Ntrials;
    
    fs = sprintf('Total error (bias & precision)'); disp(fs)
    fs = sprintf('Radial:  %.3f px\t Centroid %.3f px\t Gaussian NLLS %.3f px\t ', ...
        toterror(1,j), toterror(2,j), toterror(3,j)); disp(fs)

end

% Correlations
if length(SNrarray)==1; 
    % Correlation of residuals for radial symmetry and Gaussian methods
    % (To be plotted later)
    % Calculate...
    dx = x0_r-(N+1)/2-x0'; dy = y0_r-(N+1)/2-y0'; err_r = sqrt(dx.*dx + dy.*dy);
    dx = x0_c-(N+1)/2-x0'; dy = y0_c-(N+1)/2-y0'; err_c = sqrt(dx.*dx + dy.*dy);
%    dx = x0_g-(N+1)/2-x0'; dy = y0_g-(N+1)/2-y0'; err_g = sqrt(dx.*dx + dy.*dy);
    disp('To make plots, run make_singleSNR_plots.m')
else
    disp('To make plots, run make_multiSNr_plots.m')
end

% Save output 
x0c = x0 + (N+1)/2;
y0c = y0 + (N+1)/2;
save('test_image','z','x0c','y0c');

% Clear all the test images (z), and other things
clear j k s z

save(outfilename);



%% --------------

    function [bias sigbias sigma toterror] = calcbias(x0, y0, x0_fit, y0_fit, N)
        
        % The error of the fits
        
        % total error in x, y
        % Note that the fitting functions return the center relative to the
        % upper left pixel, so that pixel (N+1)/2 is the image center.
        dx = x0_fit-x0'-(N+1)/2;
        dy = y0_fit-y0'-(N+1)/2;
        
        % Bias
        [Ax, sigA, biasx, sigbiasx] = fitline(x0, dx, false);
        [Ay, sigA, biasy, sigbiasy] = fitline(y0, dy, false);
        bias = -0.5*(biasx+biasy);  % average in x, y
        sigbias = 0.5*(sigbiasx+sigbiasy) + 0.5*abs(sigbiasx-sigbiasy);
        % estimate uncertainty as sum of fitting and uncertainty and
        % different in directions; not perfect
        
        % Precision
        std_x = std(dx - Ax - biasx.*x0');
        std_y = std(dy - Ay - biasy.*y0');
        sigma = sqrt(std_x*std_x + std_y*std_y); % precision of the tracking error
        
        % total error (bias and precision)
        toterror = sqrt(std(dx)*std(dx) +std(dy)*std(dy));
        
    end

end

