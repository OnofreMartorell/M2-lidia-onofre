clearvars;
clc
addpath('../Block1');
addpath(genpath('.'));
%%
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
param.maxIter = 1000;
param.delta_t = 1/4;

figure(1)
imshow(I);
title('Before')


Iinp = G7_inpainting_color( I, mask, param );


figure(2)
imshow(Iinp)
title('After'); 

%% Challenge image. (We have lost 99% of information)

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
param.maxIter = 2000;
param.delta_t = 1/4;

% figure(1)
% imshow(I);
% title('Before')

Iinp = G7_inpainting_color( I, mask, param );
 
% figure(2)
imshow(Iinp)
title('After');

%% Barbara image
clearvars;

%Read the image
I = double(imread('Barbara.png'));
[ni, nj, ~] = size(I);


I = I - min(I(:));
I = I / max(I(:));

%We want to inpaint those areas in which mask == 1 (green part of the image)
I_ch1 = I(:,:,1);
I_ch2 = I(:,:,2);
I_ch3 = I(:,:,3);

%TO COMPLETE 1

mask = zeros(size(I)) ; %mask_img(i,j) == 1 means we have lost information in that pixel
                                      %mask(i,j) == 0 means we have information in that pixel
for i=1:ni
    for j=1:nj
        if(I_ch1(i,j) == 0 && I_ch2(i,j) == 1 && I_ch3(i,j) == 0)
            mask(i, j, 1) = 1;
            mask(i, j, 2) = 1;
            mask(i, j, 3) = 1;
        end
    end
end

%parameters
param.hi = 1 / (ni-1);
param.hj = 1 / (nj-1);

param.r = 10^-3;
param.lambda = 6;
param.maxIter = 5000;
param.delta_t = 1/4;


figure(1)
imshow(I);
title('Before')


Iinp = G7_inpainting_color( I, mask, param );
 
figure(2)
imshow(Iinp)
title('After');
imwrite(Iinp, './Inpaint/Restored/Barbara_restored.png');
%%

name = {'99', '95', '90', '85', '80'};
i = 5;

I = double(imread(strcat('./Inpaint/Baboon/Baboon_', name{i}, '.png')));


I = I - min(I(:));
I = I / max(I(:));

mask = imread(strcat('./Inpaint/Baboon/Baboon_', name{i}, '_mask.png'));
mask = mask > 128;


%Parameters
param.r = 10^-4;
param.lambda = 10;
param.maxIter = 1000;
param.delta_t = 1/4;


% figure(1)
% imshow(I);
% title('Before')


Iinp = G7_inpainting_color( I, mask, param );
 
imwrite(Iinp, strcat('./Inpaint/Baboon_restore/Baboon_', name{i}, '.png'))

% figure(2)
% imshow(Iinp)
% title('After');

%% Boat image
clearvars;

%Read the image
I = double(imread('Boats_modified.png'));
[ni, nj, nC] = size(I);


I = I - min(I(:));
I = I / max(I(:));

%We want to inpaint those areas in which mask == 1 (blue part of the image)

mask = imread('mask_Boats.png'); %mask_img(i,j) == 1 means we have lost information in that pixel
mask = mask > 128;                           %mask(i,j) == 0 means we have information in that pixel


% Parameters
param.r = 10^-7;
param.lambda = 7;
param.maxIter = 2000;
param.delta_t = 1/4;


% figure(1)
% imshow(I);
% title('Before')


Iinp = G7_inpainting_color( I, mask, param );
 
figure(2)
imshow(Iinp)
title('After');

imwrite(Iinp, './Inpaint/Restored/Boats_restored.png');

