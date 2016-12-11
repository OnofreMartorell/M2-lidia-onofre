name = {'99', '95', '90', '85', '80'};
i = 3;

I = double(imread(strcat('./Inpaint/Baboon/Baboon_', name{i}, '.png')));


I = I - min(I(:));
I = I / max(I(:));

mask = imread(strcat('./Inpaint/Baboon/Baboon_', name{i}, '_mask.png'));
mask = mask > 128;


%Parameters
param.r = 10^-3;
param.lambda = 10;
param.maxIter = 1000;
param.delta_t = 1/8;
param.iters_inverse = 10;



Iinp_list = inverse_scale( I, mask, param );
 
% for j = 1:param.iters_inverse
%     imwrite(uint8(Iinp_list(:, :, :, j)), strcat('./Inverse_scale_', name{i}, '_Baboon/Baboon_iter_', num2str(j), '.png'))
% end    
