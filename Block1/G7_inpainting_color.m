function [ Iinp ] = G7_inpainting_color( Im, mask )
lambda = 1/2;
delta_t = 1/4;
%Im = f
dim = size(Im);

p_n = zeros(3, 2, dim(1), dim(2));
while 1
    [ div ] = divergence( p_n(1, 1,: ,: ), p_n(1, 2,: ,: ) );
    [ grad_v_x, grad_v_y ] = gradient( div - (Im(1, :, :)/lambda) );
    grd_div_1 = [ grad_v_x, grad_v_y ];
    
    [ div ] = divergence( p_n(2, 1,: ,: ), p_n(2, 2,: ,: ) );
    [ grad_v_x, grad_v_y ] = gradient( div - (Im(2, :, :)/lambda) );
    grd_div_2 = [ grad_v_x, grad_v_y ];
    
    [ div ] = divergence( p_n(3, 1,: ,: ), p_n(3, 2,: ,: ) );
    [ grad_v_x, grad_v_y ] = gradient( div - (Im(3, :, :)/lambda) );
    grd_div_3 = [ grad_v_x, grad_v_y ];
    for i = 1:3
        p_n1 = (p_n(i, :, :, :) + delta_t )/( 1);
    end   
    
end    
end

