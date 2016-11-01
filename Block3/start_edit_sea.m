addpath('Poisson_editing_sea');



%%
% Family = double(rgb2gray(imread('family_crop.jpg')));
% Mask_family = Family < 251;
% imshow(Mask_family)
% imwrite(Mask_family, 'Poisson_editing_sea/Mask_family_bad.png');
%%
% Dog = double(rgb2gray(imread('dog_crop.jpg')));
% Mask_dog = Dog < 251;
% imshow(Mask_dog)
% imwrite(Mask_dog, 'Poisson_editing_sea/Mask_dog_bad.png');
%%
% Family = double(rgb2gray(imread('family.jpg')));
% Mask_family = rgb2gray(imread('Mask_family.png'));
% Mask_family = Mask_family > 128;
% crop = Mask_family.*Family;
% imshow(uint8(crop))

%%
% Dog = double(rgb2gray(imread('dog.jpg')));
% Mask_dog = rgb2gray(imread('Mask_dog.png'));
% Mask_dog = Mask_dog > 128;
% crop = Mask_dog.*Dog;
% imshow(uint8(crop))

%%
% Sea = double(rgb2gray(imread('the_sea.jpg')));
% Mask_family = rgb2gray(imread('Mask_family_sea.png'));
% Mask_family = Mask_family > 128;
% crop_f = Mask_family.*Sea;
% 
% Mask_dog = rgb2gray(imread('Mask_dog_sea.png'));
% Mask_dog = Mask_dog > 128;
% crop_d = Mask_dog.*Sea;
% imshow(uint8(crop_f + crop_d))
%%

clearvars;
dst = double(imread('the_sea.jpg'));
src1 = double(imread('family.jpg')); 
src2 = double(imread('dog.jpg')); 
[ni, nj, nChannels] = size(dst);

param.hi = 1;
param.hj = 1;


%masks to exchange: family swimming
mask_src = logical(rgb2gray(imread('Mask_family.png')));
mask_dst = logical(rgb2gray(imread('Mask_family_sea.png')));



for nC = 1:nChannels
    
    drivingGrad_i = sol_DiFwd(src1(:,:,nC), param.hi);
    drivingGrad_j = sol_DjFwd(src1(:,:,nC), param.hj);

    driving_on_src = sol_DiBwd(drivingGrad_i, param.hi) + sol_DjBwd(drivingGrad_j, param.hj);
    
    driving_on_dst = zeros(size(dst(:,:,1)));   
    driving_on_dst(mask_dst(:)) = driving_on_src(mask_src(:));
    
    param.driving = driving_on_dst;

    dst1(:,:,nC) = G7_Poisson_Equation_Axb(dst(:,:,nC), mask_dst,  param);
end


%Masks to exchange: Dog swimming
mask_src=logical(rgb2gray(imread('Mask_dog.png')));
mask_dst=logical(rgb2gray(imread('Mask_dog_sea.png')));
for nC = 1: nChannels
    
    %TO DO: COMPLETE the ??
    drivingGrad_i = sol_DiFwd(src2(:,:,nC), param.hi);
    drivingGrad_j = sol_DjFwd(src2(:,:,nC), param.hj);

    driving_on_src = sol_DiBwd(drivingGrad_i, param.hi) + sol_DjBwd(drivingGrad_j, param.hj);
    
    driving_on_dst = zeros(size(dst(:,:,1)));  
    driving_on_dst(mask_dst(:)) = driving_on_src(mask_src(:));
    
    param.driving = driving_on_dst;

    dst1(:,:,nC) = G7_Poisson_Equation_Axb(dst1(:,:,nC), mask_dst,  param);
end

imshow(dst1/256)



