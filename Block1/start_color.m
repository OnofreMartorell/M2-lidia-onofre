clearvars;
clc
%There are several black and white images to test:
%  image1_toRestore.jpg
%  image2_toRestore.jpg
%  image3_toRestore.jpg
%  image4_toRestore.jpg
%  image5_toRestore.jpg


name= 'image5';

I = double(imread([ name '_toRestore.jpg']));

%Normalize values into [0,1]
I = I - min(I(:));
I = I/max(I(:));

%Load the mask
mask_img = double(imread([name '_mask.jpg']));
%We want to inpaint those areas in which mask == 1
mask = mask_img >128; %mask(i,j) == 1 means we have lost information in that pixel
                      %mask(i,j) == 0 means we have information in that
                      %pixel
                                                                   
%%Parameters 
param.r = 10^-4;
param.lambda = 1/2;
param.maxIter = 100000000;

% figure(1)
% imshow(I);
% title('Before')


% Iinp = G7_inpainting_color( I, mask, param );


% figure(2)
% imshow(Iinp)
% title('After'); 

%% Challenge image. (We have lost 99% of information)
clearvars
I = double(imread('image6_toRestore.tif'));
%Normalize values into [0,1]
I = I/256;


mask_img = double(imread('image6_mask.tif'));
mask = mask_img >128; %mask(i,j) == 1 means we have lost information in that pixel
                      %mask(i,j) == 0 means we have information in that
                      %pixel

%%Parameters 
param.r = 10^-5;
param.lambda = 4;
param.maxIter = 10000;


% figure(1)
% imshow(I);
% title('Before')

Iinp = G7_inpainting_color( I, mask, param );
 
figure(2)
imshow(Iinp)
title('After');