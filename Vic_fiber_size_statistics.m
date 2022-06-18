clear; close all; clc;
[filename, pathname]=uigetfile({'G:\PhD, PMMH, ESPCI\Experimental Data (EXTRACTED)\20220614-SU8_Fibers_Size\*.tif'}, 'Choose a file to be processed');  % input file
% path for the result files
pathout = uigetdir('G:\PhD, PMMH, ESPCI\Processing\20220614-SU8_Fibers_Size\results\', 'Choose the saving folder');
[status,msg,msgID] = mkdir(pathout);

tifpath = [pathname,filesep, filename];
% gets image stack information
InfoImage=imfinfo(tifpath);
% bit depth
bitimg = InfoImage.BitDepth;
% number of images in the multi-tiff file
imtot=length(InfoImage);

for ii = 1:imtot
    Img = imread(tifpath, ii);
    fiber_img = fibermetric(Img,20,'StructureSensitivity', 0.05*diff(getrangefromclass(Img)));figure; imshow(fiber_img)

    % apply gaussian blur
    blur_img = gaussian_blur(fiber_img,10);figure; imshow(blur_img)
%     blur_img = blur_img(lzero+1:end-lzero,lzero+1:end-lzero);imshow(blur_img)

    % binarization and suppression of stray pixels
    BI = imbinarize(blur_img,'adaptive','ForegroundPolarity','bright','Sensitivity',0.1);figure; imshow(BI)
    BI = bwareafilt(BI,1); % extract object based on area, where FilNum is the expected # of filaments
    % BI = imfill(BI,'holes'); % fill the holes in the binary image of the filament (works better than bwmorph)

    % smooth out the jags/irregularities along the boundary of the objects:
    bnd = bwboundaries(BI,'noholes'); % find boundary coordinates of all objects in the FOV
    % smooth the edges of each objects in the FOV, then recontruct a binary image for each objects
    for i = 1 : FilNum
        clear xc yc
        xc = smooth(bnd{i}(:,1),1);
        yc = smooth(bnd{i}(:,2),1);
        s{i} = transpose(poly2mask(xc,yc,size(imgn,2)-2*lzero,size(imgn,1)-2*lzero));
    end
    % sum all objects, with smoothed edges, to recontruct the original binary image
    BI = logical(sum(cat(3,s{:}),3));
    % skeletonization
    skel = bwskel(BI,'MinBranchLength',MinBranchLength);

end
