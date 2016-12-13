clearvars;
addpath(genpath('../Block3'));
addpath(genpath('.'));
% dst = double(imread('lena.png'));
% src = double(imread('girl.png')); % flipped girl, because of the eyes

% dst = double(imread('the_sea.jpg'));
% src = double(imread('family.jpg')); 

dst = double(imread('madera.jpg'));
src = double(imread('firma.jpg')); 
[ni,nj, nChannels]=size(dst);

param.hi = 1;
param.hj = 1;


%Lena example
% masks_source = {'mask_src_eyes.png' 'mask_src_mouth.png'};
% masks_dest = {'mask_dst_eyes.png' 'mask_dst_mouth.png'};

masks_source = {'mask_firma.png'};
masks_dest = {'mask_madera.png'};
%masks to exchange: Eyes and Mouth
dst1 = dst;
for i = 1:length(masks_source)

    mask_src = imread(masks_source{i});
    mask_dst = imread(masks_dest{i});
    [~, ~, channels] = size(mask_src);
    if channels == 3
        mask_src = logical(rgb2gray(mask_src));
        mask_dst = logical(rgb2gray(mask_dst));
    else
        mask_src = logical(mask_src);
        mask_dst = logical(mask_dst);
    end    
    for nC = 1: nChannels
        %g: origin (source) image
        
        %Compute gradient of g and its norm
        %We want that the gradients are in the same position
        drivingGrad_i_g = sol_DiFwd(src(:,:,nC), param.hi);
        Gradient_i_g = zeros(size(dst(:,:,1)));
        Gradient_i_g(mask_dst(:)) = drivingGrad_i_g(mask_src(:));
        
        drivingGrad_j_g = sol_DjFwd(src(:,:,nC), param.hj);
        Gradient_j_g = zeros(size(dst(:,:,1)));
        Gradient_j_g(mask_dst(:)) = drivingGrad_j_g(mask_src(:));
        
        
        norm_g = Gradient_i_g.^2 + Gradient_j_g.^2;
        
        
        %f* destination image
        
        %Compute gradient of f* and its norm
        Gradient_i_f = sol_DiFwd(dst(:,:,nC), param.hi);
        Gradient_j_f = sol_DjFwd(dst(:,:,nC), param.hj);
        norm_f = Gradient_i_f.^2 + Gradient_j_f.^2;
        
        %Compute the laplacian of f*
%         driving_on_src_f = sol_DiBwd(drivingGrad_i_f, param.hi) +...
%             sol_DjBwd(drivingGrad_j_f, param.hj);
        
        greatest = norm_f > norm_g;
        
        drivingGrad_i_dst = greatest.*Gradient_i_f + (1 - greatest).*Gradient_i_g;
        drivingGrad_j_dst = greatest.*Gradient_j_f + (1 - greatest).*Gradient_j_g;
        
        driving_on_dst = sol_DiBwd(drivingGrad_i_dst, param.hi) +...
             sol_DjBwd(drivingGrad_j_dst, param.hj);
            
        param.driving = driving_on_dst;
        
        dst1(:,:,nC) = G7_Poisson_Equation_Axb(dst1(:,:,nC), mask_dst,  param);
    end
end

figure, imshow(dst1/256)

imwrite(dst1/256, strcat('Mixing/Firma.png'))

