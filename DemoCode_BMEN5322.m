% This is the demo code

clear; close all; clc;

%% Specify data file directory
% Establish filename convention used in image series

fileFolder = fullfile(pwd, 'demo');
files = dir(fullfile(fileFolder, '*dcm'));
fileNames = {files.name};

%% Extract file header (metadata, from DICOM stack)
info = dicominfo(fullfile(fileFolder, fileNames{1}));
display(info);
% Extract size info from metadata
voxel_size = [info.PixelSpacing; info.SliceThickness]';
display(voxel_size);    
% Get file information.
numImages = length(fileNames);
display(numImages);
%% Read slice images; populate 3-D matrix
%Create array
mri = zeros(info.Rows, info.Columns, numImages, class(info.BitsAllocated));
 
% for i = length(fileNames):-1:1
for i = 1:length(fileNames)
   fname = fullfile(fileFolder, fileNames{i}); 
   mri(:,:,i) = uint16(dicomread(fname));  
end


%% Creat Montage
figure;
montage(reshape(uint16(mri), [size(mri,1), size(mri,2),1, size(mri,3)]),'DisplayRange',[]); %  montage(reshape(SubV,[size(SubV,1), size(SubV,2), 1, size(SubV,3)]),'DisplayRange',[]);
set(gca, 'clim', [0, 500]); % got this value by calling i
drawnow;
shg;

D = zeros(info.Rows, info.Columns, numImages, class(info.BitsAllocated));
D(:,:,1:32) = mri(:,:,17:end);
D(:,:,33:end) = mri(:,:,1:16);
figure;
montage(reshape(uint16(D), [size(D,1), size(D,2),1, size(D,3)]),'DisplayRange',[]); %  montage(reshape(SubV,[size(SubV,1), size(SubV,2), 1, size(SubV,3)]),'DisplayRange',[]);
set(gca, 'clim', [0, 500]); % got this value by calling i
drawnow;
shg;

%% Explore image data using Image Viewer GUI toll
im = D(:,:,24); % Pick middle slice for viewing
max_level = double(max(im(:)));
imt = imtool(im, [0, max_level]);
imtool close all


%% Threshold 
lb = 50;  % lower threshold (ignore CSF & air)
ub = 300; % upper threshold (ignore subcutaneous fat);

% Use morphological operations and blob analysis tool
% Convert image to B/W and remove small blobs
mriAdjust = D;
mriAdjust(mriAdjust <= lb) = 0;
mriAdjust(mriAdjust >= ub) = 0;
mriAdjust(:,:,1:21) = 0;
mriAdjust(:,:,43:48) = 0;
bw = logical(mriAdjust); % convert to binary

%% Let's clean that up a bit
nhood = ones([5 5 1]);
bw = imopen(bw,nhood);
figure;
subplot(1,2,1);imshow(D(:,:,24),[0, max_level]);
subplot(1,2,2);imshow(bw(:,:,24));



%% lesion area segmentation
% Identify blobs and display them for a single image
L        = bwlabeln(bw);
stats    = regionprops(L, 'Area','Centroid');
% Determin the largest blob and eliminate all others
A        = [stats.Area];
biggest  = find(A == max(A));
mriAdjust(L ~= biggest) = 0;

imA      = imadjust(mriAdjust(:,:,24));
figure;
imshow(imA);
figure;
montage(reshape(uint16(mriAdjust), [size(mriAdjust,1), size(mriAdjust,2),1, size(mriAdjust,3)]),'DisplayRange',[]); %  montage(reshape(SubV,[size(SubV,1), size(SubV,2), 1, size(SubV,3)]),'DisplayRange',[]);
drawnow;
shg;


%% Partition brain mass
% Separate the brain mass into 2 categories
thresh_tool(uint16(mriAdjust(:,:,24)), 'gray'); % Need to download this tool - set threshold level

level = 200;

mriBrainPartition = uint8(zeros(size(mriAdjust)));
mriBrainPartition(mriAdjust>lb & mriAdjust<level) = 2;
mriBrainPartition(mriAdjust>=level) = 20;
figure;
imshow(mriBrainPartition(:,:,24), [0 0 0; 0.1 0.1 0.1; 0.25 0.25 0.25;0.25 0.25 0.8]);
figure;
montage(mriBrainPartition,'Size', [7 7],'DisplayRange', [0 20]);

%% Contour of the nth slice?
cm = brighten(jet(60),-0.5);
figure('Colormap', cm)
contourslice(mriAdjust, [], [], 24)
axis ij tight
daspect([1,1,1])













