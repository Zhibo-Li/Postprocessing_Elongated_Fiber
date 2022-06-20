clear; close all; clc;
[filename, pathname]=uigetfile({'C:\Users\LaVision\Desktop\20220614-SU8_Fibers_Size\*.tif'}, 'Choose a file to be processed');  % input file
% path for the result files
% pathout = uigetdir('C:\Users\LaVision\Desktop\20220614-SU8_Fibers_Size\', 'Choose the saving folder');
% [status,msg,msgID] = mkdir(pathout);

tifpath = [pathname,filesep, filename];
% gets image stack information
InfoImage=imfinfo(tifpath);
% bit depth
bitimg = InfoImage.BitDepth;
% number of images in the multi-tiff file
imtot=length(InfoImage);


% define the parameters of func:Vic_get_skeleton.
thickness = 20;
structsens = 0.1;
lnoise = 10;
lobject = false;
threshold = 0;
bina_sensitivity = 0.1;
MinBranchLength = 300;
FilNum = 3;

for img_i = 5%1:imtot
    Img = imread(tifpath, img_i);
    IMG = im2double(Img) - 0.9;
    imshow(Img); figure; imshow(uint16(rescale(IMG, 0, 2^bitimg-1)));

    [fiber_img,blur_img,BI,skel_full,L_full,cntrds] = Vic_get_skeleton(Img, ...
        thickness,structsens,lnoise,lobject,threshold,bina_sensitivity,MinBranchLength,FilNum);

end


function [fiber_img,blur_img,BI,skel_full,L_full,cntrds] = Vic_get_skeleton(ImageIN,...
    thickness,structsens,lnoise,lobject,threshold,bina_sensitivity,MinBranchLength,FilNum)
%
% NAME:
%               Vic_get_skeleton
% PURPOSE:
%               Find the skeletons of the elongated shapes.
% CATEGORY:
%               Image Processing
% CALLING SEQUENCE:
%               [fiber_img,blur_img,BI,skel_full,L_full,cntrds] = Vic_get_skeleton(
%               ImageIN,thickness,structsens,lnoise,lobject,threshold,bina_sensitivity,
%               MinBranchLength,FilNum)
%
% INPUTS:
%               ImageIN:    The two-dimensional array to be processed.
%               thickness and structsens:   SEE func:fibermetric.
%                       ! The value of structsens indicates the percentage
%                       of the diff(getrangefromclass(I))
%               lnoise, lobject and threshold:  SEE func:gaussian_blur.
%                       lobject and threshold are optional, setting to:
%                       lobject = false. (by default)
%                       threshold = 0 (by default).
%               bina_sensitivity:   SEE func:imbinarize.
%               MinBranchLength:    SEE func:bwskel.
%               FilNum: SEE func:bwareafilt.
%                       The number of the filaments you expected.
%
% OUTPUTS:
%               fiber_img, blur_img, BI, skel_full, L_full, cntrds.
% NOTES:

% MODIFICATION HISTORY:
%               Original from Francesco BONACCI, ESPCI Paris, 2020.
%
%               Written by Zhibo LI, ESPCI Paris, 2022-06.

% Enhance elongated or tubular structures in images
fiber_img = fibermetric(ImageIN, thickness, 'StructureSensitivity', ...
    structsens*diff(getrangefromclass(ImageIN)));
%figure;imshow(fiber_img)

% apply gaussian blur
blur_img = gaussian_blur(fiber_img, lnoise, lobject, threshold); %figure;imshow(blur_img)

% binarization and suppression of stray pixels
BI = imbinarize(blur_img, 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', bina_sensitivity);
BI = bwareafilt(BI, FilNum); % extract object based on area, where FilNum is the expected # of filaments
BI = imfill(BI, 'holes'); % fill the holes in the binary image of the filament (works better than bwmorph) 
 %figure;imshow(BI)
bnd = bwboundaries(BI, 'noholes'); % find boundary coordinates of all objects in the FOV
% smooth the edges of each objects in the FOV, then recontruct a binary image for each objects
for filament_i = 1 : FilNum
    clear xc yc
    xc = smooth(bnd{filament_i}(:,1),1);
    yc = smooth(bnd{filament_i}(:,2),1);
    s_full{filament_i} = transpose(poly2mask(xc,yc,size(ImageIN,2),size(ImageIN,1)));
end
BI_full = logical(sum(cat(3,s_full{:}),3));
skel_full = bwskel(BI_full,'MinBranchLength', MinBranchLength); 
figure;imshow(labeloverlay(ImageIN, skel_full, 'Transparency', 0))
% find the connected part of the skeletonization image
CC_full = bwconncomp(skel_full);
% label each filament with a different number
L_full = labelmatrix(CC_full);
% find the centroids
cntrds = regionprops(L_full ,'centroid');

end