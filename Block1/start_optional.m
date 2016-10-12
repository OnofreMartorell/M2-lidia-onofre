%% Woman Image
clearvars;

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
title('After');

%% Boat Image
clearvars;

%Read the image
I = double(imread('Boats24.tiff'));

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
        if(I_ch1(i,j) == 0 && I_ch2(i,j) == 0 && I_ch3(i,j) == 1)
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