addpath('Poisson_editing_motos');



%%
% Thunderbolt = double(rgb2gray(imread('Thunderbolt_crop.png')));
% Mask_thunderbolt = Thunderbolt < 251;
% imshow(Mask_thunderbolt)
% imwrite(Mask_thunderbolt, 'Poisson_editing_motos/Mask_thunderbolt_bad.png');

%%
% Thunderbolt = double(rgb2gray(imread('thunderbolt.png')));
% Mask_thunderbolt = rgb2gray(imread('Mask_thunderbolt.png'));
% Mask_thunderbolt = Mask_thunderbolt > 128;
% crop = Mask_thunderbolt.*Thunderbolt;
% imshow(uint8(crop))


%%
% Sea = double(rgb2gray(imread('Mar_motos.jpg')));
% 
% Mask_thunderbolt01 = rgb2gray(imread('Mask_motos_01.jpg'));
% Mask_thunderbolt01 = Mask_thunderbolt01 > 128;
% crop_01 = Mask_thunderbolt01.*Sea;
% 
% Mask_thunderbolt02 = rgb2gray(imread('Mask_motos_02.jpg'));
% Mask_thunderbolt02 = Mask_thunderbolt02 > 128;
% crop_02 = Mask_thunderbolt02.*Sea;
% 
% Mask_thunderbolt03 = rgb2gray(imread('Mask_motos_03.jpg'));
% Mask_thunderbolt03 = Mask_thunderbolt03 > 128;
% crop_03 = Mask_thunderbolt03.*Sea;
% 
% Mask_thunderbolt04 = rgb2gray(imread('Mask_motos_04.jpg'));
% Mask_thunderbolt04 = Mask_thunderbolt04 > 128;
% crop_04 = Mask_thunderbolt04.*Sea;
% 
% imshow(uint8(crop_01 + crop_02 + crop_03 + crop_04))
%%

clearvars;
dst = double(imread('Mar_motos.jpg'));
src = double(imread('Thunderbolt.png')); 

[ni, nj, nChannels] = size(dst);

param.hi = 1;
param.hj = 1;



mask_src = logical(rgb2gray(imread('Mask_thunderbolt.png')));
mask_dst = logical(rgb2gray(imread('Mask_motos_01.png')));



for nC = 1:nChannels
    
    drivingGrad_i = sol_DiFwd(src(:,:,nC), param.hi);
    drivingGrad_j = sol_DjFwd(src(:,:,nC), param.hj);

    driving_on_src = sol_DiBwd(drivingGrad_i, param.hi) + sol_DjBwd(drivingGrad_j, param.hj);
    
    driving_on_dst = zeros(size(dst(:,:,1)));   
    driving_on_dst(mask_dst(:)) = driving_on_src(mask_src(:));
    
    param.driving = driving_on_dst;

    dst1(:,:,nC) = G7_Poisson_Equation_Axb(dst(:,:,nC), mask_dst,  param);
end

imshow(dst1/256)



