% sif2tif.m
% 
% Transform Andor *.sif file into *.tif image stack, currently only
% videos with certain resolution and length is accepted.
% 
%
% INPUT:
%       filename: the name of the sif file
%       framenumber: total number of frames in the sif file
%       Width: the width of the image
%       Height: the height of the image
%       resolution: an 1*2 vector, often [width, height]
%       NumFrames: specific index of frame to transform
%       
% OUTPUT:
%       imageData: image of size defined by the input resolution
%       
% 
%
% 
% Tzu-Chi Yen
% last modified Sep 2014
% began March 2014
% 
% Modified from Uli Kleinger's version in MATLAB CENTRAL: #11224
% Originally by Marcel Leutenegger November 2006




function imageData = sif2tif(filename, framenumber, Width, Height, resolution, NumFrames)

file = filename; 
currentFrameNumber = framenumber;  
f=fopen(file,'r');
if f < 0 
error('Could not open the file.'); 
end 
if ~isequal(fgetl(f),'Andor Technology Multi-Channel File') 
fclose(f); 
error('Not an Andor SIF image file.'); 
end 
skipLines(f,1); 
[imageData]=readSection(f, currentFrameNumber, Width, Height, resolution, NumFrames); 
fclose(f);
if currentFrameNumber > 0 
pic_out = imageData; 
%pic_plot_uint8 = uint8(pic_out/max(pic_out(:))*255); 
pic_uint8 = uint8(round(pic_out/(2^16)*65536)); 
%asdasd 
%size(pic_uint8) 

end

%figure 
%imagesc(pic_uint8_rgb)

%Read a file section. 
% 
% f File handle 
% info Section data 
% next Flags if another section is available 
% 
function [imageData]=readSection(f, currentFrameNumber, Width, Height, resolution, NumFrames)

settingsline = fgetl(f); 
settings_cell = regexp(settingsline, '(\S)+\s', 'match'); 

%% Werte aus dem Sif file 
% info.date = datestr(timestamp_total_matlab_corrected); 


%% Daraus berechnete Werte 

%% In die Ausgabe uebernehmen 

%% Weitere Daten aus dem Sif file  

current_line_text = ''; 
while isempty(strfind(current_line_text, 'Pixel number')) 
current_line_text = fgetl(f); 
end 
skipLines(f,1) 
imageArealine = fgetl(f); 
frameArealine = fgetl(f); 
imageAreaData=sscanf(imageArealine,'Pixel number%d %d %d %d %d %d %d %d %d'); 
frameAreaData=sscanf(frameArealine,'%d %d %d %d %d %d %d'); 
pixelPerImage = imageAreaData(9); 
pixelInVideo = imageAreaData(8); 
leftPixel = frameAreaData(2); 
topPixel = frameAreaData(3); 
rightPixel = frameAreaData(4); 
bottomPixel = frameAreaData(5); 



%Width=128;
%Height=128;
%resolution = [128 128];
%NumFrames = 25000;

if currentFrameNumber >0 && currentFrameNumber <= NumFrames 
% old trail-and-error traces
%num_bytes_to_skip_for_curr_frame = 4 * ((currentFrameNumber-1) * prod(resolution)); 
%num_bytes_to_skip_for_curr_frame = (4*(9+currentFrameNumber)-1)*prod(resolution)+8*resolution(1);
%num_bytes_to_skip_for_curr_frame = (4*(9+currentFrameNumber)-1)*prod(resolution)+8*resolution(1)+1250*resolution(1);

%=============================================================
% Change the parameter settings if you are transforming *.sif files 
% that are not 128*128
%=============================================================


% size: 128; NumFrames: any
num_bytes_to_skip_for_curr_frame = (4*(9+currentFrameNumber)-1)*prod(resolution)+8*resolution(1)+(NumFrames-20000)*(1/4)*resolution(1);

% size: 96; NumFrames: 50000
%num_bytes_to_skip_for_curr_frame = (4*(23+currentFrameNumber)-1)*prod(resolution)+8*resolution(1)+(NumFrames-20000)*(1/4)*resolution(1)+928*4;

% size: 64; NumFrames: 10000
%num_bytes_to_skip_for_curr_frame = (4*(28+currentFrameNumber)-1)*prod(resolution)+8*resolution(1)+(NumFrames-20000)*(1/4)*resolution(1)+2112*4;

% size: 32; NumFrames: 200000
%num_bytes_to_skip_for_curr_frame = (4*(28+currentFrameNumber)-1)*prod(resolution)+8*resolution(1)+(NumFrames-20000)*(1/4)*resolution(1)+2112*4+1180*32^2*4+64*4;


fseek(f, num_bytes_to_skip_for_curr_frame, 'cof');
fseek(f, 2, 'cof'); %keine ahnung warum 2 byte offset.. 
if prod(resolution) ~= pixelPerImage || pixelPerImage*NumFrames ~= pixelInVideo 
fclose(f) 
error('Inconsistent image header.'); 
end

imageData=(reshape(fread(f,prod(resolution),'single=>single'), Width, Height ))'; 
end

%Read a character string. 
% 
% f File handle 
% o String 
% 
function o=readString(f) 
n=fscanf(f,'%d',1); 
if isempty(n) || n < 0 || isequal(fgetl(f),-1) 
fclose(f); 
error('Inconsistent string.'); 
end 
o=fread(f,[1 n],'uint8=>char');

%Read a line. 
% 
% f File handle 
% o Read line 
% 
function o=readLine(f) 
o=fgetl(f); 
if isequal(o,-1) 
fclose(f); 
error('Inconsistent image header.'); 
end 
o=deblank(o);

%Skip bytes. 
% 
% f File handle 
% N Number of bytes to skip 
% 
function skipBytes(f,N) 
[~,n]=fread(f,N,'uint8'); 
if n < N 
fclose(f); 
error('Inconsistent image header.'); 
end

%Skip lines. 
% 
% f File handle 
% N Number of lines to skip 
% 
function skipLines(f,N) 
for n=1:N 
if isequal(fgetl(f),-1) 
fclose(f); 
error('Inconsistent image header.'); 
end 
end