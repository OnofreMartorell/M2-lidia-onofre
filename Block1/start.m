%Example script: You should replace the beginning of each function ('sol')
%with the name of your group. i.e. if your gropu name is 'G8' you should
%call :
% G8_DualTV_Inpainting_GD(I, mask, paramInp, paramROF)

clearvars;
%There are several black and white images to test:
%  image1_toRestore.jpg
%  image2_toRestore.jpg
%  image3_toRestore.jpg
%  image4_toRestore.jpg
%  image5_toRestore.jpg


name= 'image5';

I = double(imread([ name '_toRestore.jpg']));
%I=I(1:10,1:10);

%Number of pixels for each dimension, and number of channles
[ni, nj, nC] = size(I);

if nC == 3
    I = mean(I,3); %Convert to b/w. If you load a color image you should comment this line
end

%Normalize values into [0,1]
I=I-min(I(:));
I=I/max(I(:));

%Load the mask
mask_img = double(imread([name '_mask.jpg']));
%mask_img =mask_img(1:10,1:10);
[~, ~, nC] = size(mask_img);
if nC == 3
    mask_img = mask_img(:,:,1); %Convert to b/w. If you load a color image you should comment this line
end
%We want to inpaint those areas in which mask == 1
mask = mask_img >128; %mask(i,j) == 1 means we have lost information in that pixel
                      %mask(i,j) == 0 means we have information in that
                      %pixel
                                                                    
%%%Parameters for gradient descent (you do not need for week1)
param.dt = 5*10^-7;
param.iterMax = 10^4;
param.tol = 10^-5;

%%Parameters 
param.hi = 1 / (ni-1);
param.hj = 1 / (nj-1);


figure(1)
imshow(I);
title('Before')


Iinp = G7_sol_Laplace_Equation_Axb(I, mask, param);


% Iinp=G7_sol_Alternative_Laplace_Equation_Axb(I, mask, param);



figure(2)
imshow(Iinp)
title('After');

%% Challenge image. (We have lost 99% of information)
clearvars
I=double(imread('image6_toRestore.tif'));
%Normalize values into [0,1]
I=I/256;


%Number of pixels for each dimension, and number of channels

[ni, nj, ~] = size(I);

mask_img=double(imread('image6_mask.tif'));
mask = mask_img >128; %mask(i,j) == 1 means we have lost information in that pixel
                      %mask(i,j) == 0 means we have information in that
                      %pixel

param.hi = 1 / (ni-1);
param.hj = 1 / (nj-1);

figure(1)
imshow(I);
title('Before')

Iinp(:,:,1)=G7_sol_Laplace_Equation_Axb(I(:,:,1), mask(:,:,1), param);
Iinp(:,:,2)=G7_sol_Laplace_Equation_Axb(I(:,:,2), mask(:,:,2), param);
Iinp(:,:,3)=G7_sol_Laplace_Equation_Axb(I(:,:,3), mask(:,:,3), param);

% Iinp(:,:,1)=G7_sol_Alternative_Laplace_Equation_Axb(I(:,:,1), mask(:,:,1), param);
% Iinp(:,:,2)=G7_sol_Alternative_Laplace_Equation_Axb(I(:,:,2), mask(:,:,2), param);
% Iinp(:,:,3)=G7_sol_Alternative_Laplace_Equation_Axb(I(:,:,3), mask(:,:,3), param);
 

figure(2)
imshow(Iinp)
title('After');

%% Goal Image
clearvars;

%Read the image
I = double(imread('image_to_Restore.png'));

i = I(:,:,1);
[ni, nj, ~] = size(I);


I = I - min(I(:));
I = I / max(I(:));

%We want to inpaint those areas in which mask == 1 (red part of the image)
I_ch1 = I(:,:,1);
I_ch2 = I(:,:,2);
I_ch3 = I(:,:,3);

%TO COMPLETE 1


%Create the mask
mask = zeros(ni,nj) ; 
%mask_img(i,j) == 1 means we have lost information in that pixel
                                      %mask(i,j) == 0 means we have information in that pixel
for i=1:ni
    for j=1:nj
%         if(I_ch1(i,j) == 1 && I_ch2(i,j) == 0 && I_ch3(i,j) == 0 && i < 284)
        if(I_ch1(i,j) == 1 && I_ch2(i,j) == 0 && I_ch3(i,j) == 0)
            mask(i,j) = 1;
        end
    end
end

%%%Parameters for gradient descent (you do not need for week1)
%param.dt = 5*10^-7;
%param.iterMax = 10^4;
%param.tol = 10^-5;

%parameters
param.hi = 1 / (ni-1);
param.hj = 1 / (nj-1);

% for each channel 

figure(1)
imshow(I);
title('Before')

Iinp(:,:,1)=G7_sol_Laplace_Equation_Axb(I_ch1, mask, param);
Iinp(:,:,2)=G7_sol_Laplace_Equation_Axb(I_ch2, mask, param);
Iinp(:,:,3)=G7_sol_Laplace_Equation_Axb(I_ch3, mask, param);
    
figure(2)
imshow(Iinp)
% title('After');
% imwrite(Iinp, 'Goal_inpaint_bad.png');
