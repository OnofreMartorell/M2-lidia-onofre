%% Woman Image
clearvars;
addpath('../Block1');
addpath(genpath('.'));
%Read the image
I = double(imread('Barbara24.tiff'));

[ni, nj, nC] = size(I);


I = I - min(I(:));
I = I / max(I(:));

%We want to inpaint those areas in which mask == 1 (red part of the image)
I_ch1 = I(:,:,1);
I_ch2 = I(:,:,2);
I_ch3 = I(:,:,3);

%TO COMPLETE 1

mask = zeros(ni,nj) ; %mask_img(i,j) == 1 means we have lost information in that pixel
                                      %mask(i,j) == 0 means we have information in that pixel
for i=1:ni
    for j=1:nj
        if(I_ch1(i,j) == 0 && I_ch2(i,j) == 1 && I_ch3(i,j) == 0)
            mask(i,j) = 1;
        end
    end
end

%parameters
param.hi = 1 / (ni-1);
param.hj = 1 / (nj-1);

% for each channel 

figure(1)
imshow(I(:,:,1:3))


title('Before')
for i = 1:3
    Iinp(:, :, i) = G7_sol_Laplace_Equation_Axb(I(:, :, i), mask, param);
end

    
figure(2)
imshow(Iinp)
title('After');

%% Boat Image
clearvars;

%Read the image
I = double(imread('Boats.png'));
pre_mask2 = rgb2gray(imread('pre-mask_Boats.png'));
pre_mask = pre_mask2 < 128;
mask = zeros(size(I));
mask(:, :, 1) = pre_mask;
mask(:, :, 2) = pre_mask;
mask(:, :, 3) = pre_mask;

modified_im = I;
modified_im(:, :, 1) = (1 - pre_mask).*I(:, :, 1);
modified_im(:, :, 2) = (1 - pre_mask).*I(:, :, 2);
modified_im(:, :, 3) = (1 - pre_mask).*I(:, :, 3) + 255*pre_mask;
 

imwrite(mask, './Inpaint/mask_Boats.png');
imwrite(uint8(modified_im), './Inpaint/Boats_modified.png');
[ni, nj, nC] = size(I);


I = I - min(I(:));
I = I / max(I(:));

%We want to inpaint those areas in which mask == 1 (red part of the image)
I_ch1 = I(:,:,1);
I_ch2 = I(:,:,2);
I_ch3 = I(:,:,3);




%%%Parameters for gradient descent (you do not need for week1)
%param.dt = 5*10^-7;
%param.iterMax = 10^4;
%param.tol = 10^-5;

%parameters
param.hi = 1 / (ni-1);
param.hj = 1 / (nj-1);

% for each channel 
% 
% figure(1)
% imshow(I(:, :, 1:3));
% title('Before')
% 
% 
% for i = 1:3
%     Iinp(:,:,i) = G7_sol_Laplace_Equation_Axb(I(:,:,i), mask, param);
% end
% 
%     
% figure(2)
% imshow(Iinp)
% title('After');

%%
I = double(imread('Baboon.png'));
P = [0.01 0.05 0.1 0.15 0.2];
name = {'99', '95', '90', '85', '80'};
for i = 1:length(P)
    values_clear = binornd(1, P(i), size(I));
    new_im = values_clear.*I;
    imwrite(uint8(new_im), strcat('./Inpaint/Baboon/Baboon_', name{i}, '.png'))
    imwrite(1 - values_clear, strcat('./Inpaint/Baboon/Baboon_', name{i}, '_mask.png'))
end
