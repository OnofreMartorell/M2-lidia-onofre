function [ Iinp_list ] = inverse_scale( Im, mask, param )


iters = param.iters_inverse;

Iinp_list = zeros([size(Im) iters]);


v_k1 = zeros(size(Im));

for k = 1:iters
    
    f_star = Im - v_k1;
    
    Iinp  = G7_inpainting_color( f_star, mask, param );
    Iinp_list(:, :, :, k) = Iinp;
    
    v_k = v_k1;
    v_k1 = v_k + Im - Iinp;
    
    imshow(Iinp)
    imwrite(Iinp, strcat('./Inverse_scale_90_Baboon/Baboon_iter_', num2str(k), '.png'))
end

