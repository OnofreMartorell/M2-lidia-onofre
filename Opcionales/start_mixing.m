clearvars;
addpath(genpath('../Block3'));
addpath(genpath('.'));
dst = double(imread('lena.png'));
src = double(imread('girl.png')); % flipped girl, because of the eyes

% dst = double(imread('the_sea.jpg'));
% src = double(imread('family.jpg')); 

% dst = double(imread('madera.jpg'));
% src = double(imread('firma.jpg')); 
[ni,nj, nChannels]=size(dst);

param.hi = 1;
param.hj = 1;


%Lena example
masks_source = {'mask_src_eyes.png' 'mask_src_mouth.png'};
masks_dest = {'mask_dst_eyes.png' 'mask_dst_mouth.png'};

% masks_source = {'mask_firma.png'};
% masks_dest = {'mask_madera.png'};
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
        drivingGrad_i_g = sol_DiFwd(src(:,:,nC), param.hi);
        drivingGrad_j_g = sol_DjFwd(src(:,:,nC), param.hj);
        norm_grad = drivingGrad_i_g.^2 + drivingGrad_j_g.^2;
        
        %We want that the gradients are in the same position
        
        norm_g = zeros(size(dst(:,:,1)));
        norm_g(mask_dst(:)) = norm_grad(mask_src(:));
        
        %Compute the laplacian of g
        driving_on_src = sol_DiBwd(drivingGrad_i_g, param.hi) + ...
            sol_DjBwd(drivingGrad_j_g, param.hj);
        driving_on_src_g = zeros(size(dst(:,:,1)));
        driving_on_src_g(mask_dst(:)) = driving_on_src(mask_src(:));
        
        %f* destination image
        
        %Compute gradient of f* and its norm
        drivingGrad_i_f = sol_DiFwd(dst(:,:,nC), param.hi);
        drivingGrad_j_f = sol_DjFwd(dst(:,:,nC), param.hj);
        norm_grad = drivingGrad_i_f.^2 + drivingGrad_j_f.^2;
        
        norm_f = zeros(size(dst(:,:,1)));
        norm_f(mask_dst(:)) = norm_grad(mask_dst(:));
        %Compute the laplacian of f*
        driving_on_src_f = sol_DiBwd(drivingGrad_i_f, param.hi) +...
            sol_DjBwd(drivingGrad_j_f, param.hj);
        
        greatest = norm_f > norm_g;
        
%         greatest = zeros(size(norm_f));
        
        driving = greatest.*driving_on_src_f + (1 - greatest).*driving_on_src_g;
        alpha = 0;
%         driving = alpha*driving_on_src_f + (1 - alpha)*driving_on_src_g;
        
        driving_on_dst = zeros(size(dst(:,:,1)));
        driving_on_dst(mask_dst(:)) = driving(mask_dst(:));
        
        param.driving = driving_on_dst;
        
        dst1(:,:,nC) = G7_Poisson_Equation_Axb(dst1(:,:,nC), mask_dst,  param);
    end
end

figure, imshow(dst1/256)

imwrite(dst1/256, strcat('Mixing/Lena_mixing.png'))

%%

% driving = zeros(size(greatest));
%         for j1 = 1:size(greatest, 1)
%             for j2 = 1:size(greatest, 2)
%                 if norm_f(j1, j2) > norm_g(j1, j2)
%                     driving(j1, j2) = driving_on_src_f(j1, j2);
%                 else
%                     driving(j1, j2) = driving_on_src_g(j1, j2);
%                 end    
%             end
%         end    
% %Mouth
% %masks to exchange: Mouth
% mask_src=logical(imread('mask_src_mouth.png'));
% mask_dst=logical(imread('mask_dst_mouth.png'));
% for nC = 1: nChannels
%     
%     %g: origin (source) image
%     
%     %Compute gradient of g and its norm
%     drivingGrad_i_g = sol_DiFwd(src(:,:,nC), param.hi);
%     drivingGrad_j_g = sol_DjFwd(src(:,:,nC), param.hj);
%     norm_grad = drivingGrad_i_g.^2 + drivingGrad_j_g.^2;
%     
%     %We want that the gradients are in the same position
%     
%     norm_g = zeros(size(src(:,:,1)));
%     norm_g(mask_dst(:)) = norm_grad(mask_src(:));
%     
%     %Compute the laplacian of g
%     driving_on_src = sol_DiBwd(drivingGrad_i_g, param.hi) + ...
%         sol_DjBwd(drivingGrad_j_g, param.hj);
%     driving_on_src_g = zeros(size(src(:,:,1)));
%     driving_on_src_g(mask_dst(:)) = driving_on_src(mask_src(:));
%     
%     %f* destination image
%     
%     %Compute gradient of f* and its norm
%     drivingGrad_i_f = sol_DiFwd(dst(:,:,nC), param.hi);
%     drivingGrad_j_f = sol_DjFwd(dst(:,:,nC), param.hj);
%     norm_grad = drivingGrad_i_f.^2 + drivingGrad_j_f.^2;
%     
%     norm_f = zeros(size(dst(:,:,1)));
%     norm_f(mask_dst(:)) = norm_grad(mask_dst(:));
%     %Compute the laplacian of f*
%     driving_on_src_f = sol_DiBwd(drivingGrad_i_f, param.hi) +...
%         sol_DjBwd(drivingGrad_j_f, param.hj);
%     
%     greatest = norm_f > norm_g;
%     
%     driving = greatest.*driving_on_src_f + (1 - greatest).*driving_on_src_g;
%     
%     
%     driving_on_dst = zeros(size(src(:,:,1)));
%     driving_on_dst(mask_dst(:)) = driving(mask_dst(:));
%     
%     param.driving = driving_on_dst;
%     
%     dst1(:,:,nC) = G7_Poisson_Equation_Axb(dst1(:,:,nC), mask_dst,  param);
% end



%%
% im = zeros(ni, nj, 3);
%
% for i = 1:3
%    I = dst(:, :, i);
%    ii = src(:, :, i);
%    I(mask_dst(:)) = ii(mask_src(:));
%    im(:, :, i) = I;
% end
% figure, imshow(uint8(im))