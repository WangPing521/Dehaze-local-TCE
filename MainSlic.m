tic
close all;
clc;
clear all;

% parameters related with the region division
K =16; % tunning parameter, the number of seeds for SLIC clusters
m_compactness = 55; % tunning parameter, SLIC algorithm related parameter

img = imread('45_hazy.jpg');
figure
imshow(img);
title('the inut hazy image');

img_size = size(img);   

%shift the image from RGB to LAB space
cform = makecform('srgb2lab');       
img_Lab = applycform(img, cform);    

img_sz = img_size(1)*img_size(2);
superpixel_sz = img_sz/K;
STEP = uint32(sqrt(superpixel_sz));
xstrips = uint32(img_size(2)/STEP);
ystrips = uint32(img_size(1)/STEP);
xstrips_adderr = double(img_size(2))/double(xstrips);
ystrips_adderr = double(img_size(1))/double(ystrips);
numseeds = xstrips*ystrips;

kseedsx = zeros(numseeds, 1);
kseedsy = zeros(numseeds, 1);
kseedsl = zeros(numseeds, 1);
kseedsa = zeros(numseeds, 1);
kseedsb = zeros(numseeds, 1);
n = 1;
for y = 1: ystrips
    for x = 1: xstrips 
        kseedsx(n, 1) = (double(x)-0.5)*xstrips_adderr;
        kseedsy(n, 1) = (double(y)-0.5)*ystrips_adderr;
        kseedsl(n, 1) = img_Lab(fix(kseedsy(n, 1)), fix(kseedsx(n, 1)), 1);
        kseedsa(n, 1) = img_Lab(fix(kseedsy(n, 1)), fix(kseedsx(n, 1)), 2);
        kseedsb(n, 1) = img_Lab(fix(kseedsy(n, 1)), fix(kseedsx(n, 1)), 3);
        n = n+1;
    end
end
n = 1;
%compute super-pixels
klabels = PerformSuperpixelSLIC(img,img_Lab, kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, STEP, m_compactness);
img_Contours = DrawContoursAroundSegments(img, klabels);
%merge small regions
nlabels = EnforceLabelConnectivity(img_Lab, klabels, K); 
%draw the contour of each region and dehaze
[light,Anum,img_ContoursEX] = DrawContoursAroundSegments_EX(img, nlabels,numseeds);
toc
        
        






