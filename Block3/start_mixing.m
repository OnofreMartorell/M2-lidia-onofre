clearvars;
dst = double(imread('lena.png'));
src = double(imread('girl.png')); % flipped girl, because of the eyes
[ni,nj, nChannels]=size(dst);

param.hi=1;
param.hj=1;


%masks to exchange: Eyes
mask_src=logical(imread('mask_src_eyes.png'));
mask_dst=logical(imread('mask_dst_eyes.png'));

for nC = 1: nChannels
    
    %TO DO: COMPLETE the ??
    drivingGrad_i_g = sol_DiFwd(src(:,:,nC), param.hi);
    drivingGrad_j_g = sol_DjFwd(src(:,:,nC), param.hj);
    
    driving_on_src_g = sol_DiBwd(drivingGrad_i_g, param.hi) +...
        sol_DjBwd(drivingGrad_j_g, param.hj);
    
    norm_g = drivingGrad_i_g.^2 + drivingGrad_j_g.^2;
    
    norm_g(mask_src(:)) = norm_g(mask_src(:));
    
    drivingGrad_i_f = sol_DiFwd(mask_dst.*dst(:,:,nC), param.hi);
    drivingGrad_j_f = sol_DjFwd(mask_dst.*dst(:,:,nC), param.hj);
    
    driving_on_src_f = sol_DiBwd(drivingGrad_i_f, param.hi) +...
        sol_DjBwd(drivingGrad_j_f, param.hj);
    
    norm_f = drivingGrad_i_f.^2 + drivingGrad_j_f.^2;
    
    norm_f(mask_dst(:)) = norm_f(mask_src(:));
    
    greatest = norm_f > norm_g;
    
    driving = greatest.*driving_on_src_f + (1 - greatest)*driving_on_src_g;
    
    driving_on_dst = zeros(size(src(:,:,1)));   
    driving_on_dst(mask_dst(:)) = driving(mask_src(:));
    
    param.driving = driving_on_dst;

    dst1(:,:,nC) = G7_Poisson_Equation_Axb(dst(:,:,nC), mask_dst,  param);
end



%Mouth
%masks to exchange: Mouth
mask_src=logical(imread('mask_src_mouth.png'));
mask_dst=logical(imread('mask_dst_mouth.png'));
for nC = 1: nChannels
    
    %TO DO: COMPLETE the ??
    drivingGrad_i_g = sol_DiFwd(src(:,:,nC), param.hi);
    drivingGrad_j_g = sol_DjFwd(src(:,:,nC), param.hj);
    
    driving_on_src_g = sol_DiBwd(drivingGrad_i_g, param.hi) +...
        sol_DjBwd(drivingGrad_j_g, param.hj);
    
    norm_g = drivingGrad_i_g.^2 + drivingGrad_j_g.^2;
    
    norm_g(mask_src(:)) = norm_g(mask_src(:));
    
    drivingGrad_i_f = sol_DiFwd(mask_dst.*dst(:,:,nC), param.hi);
    drivingGrad_j_f = sol_DjFwd(mask_dst.*dst(:,:,nC), param.hj);
    
    driving_on_src_f = sol_DiBwd(drivingGrad_i_f, param.hi) +...
        sol_DjBwd(drivingGrad_j_f, param.hj);
    
    norm_f = drivingGrad_i_f.^2 + drivingGrad_j_f.^2;
    
    norm_f(mask_dst(:)) = norm_f(mask_src(:));
    greatest = norm_f > norm_g;
    
    driving = greatest.*driving_on_src_f + (1 - greatest)*driving_on_src_g;
    
    driving_on_dst = zeros(size(src(:,:,1)));  
    driving_on_dst(mask_dst(:)) = driving(mask_src(:));
    
    param.driving = driving_on_dst;

    dst1(:,:,nC) = G7_Poisson_Equation_Axb(dst1(:,:,nC), mask_dst,  param);
end

imshow(dst1/256)

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