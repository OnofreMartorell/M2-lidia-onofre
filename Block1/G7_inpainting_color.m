function [ Iinp ] = G7_inpainting_color( Im, mask, param )
lambda = param.lambda;
delta_t = 1/4;
r = param.r;
maxIter = param.maxIter;


dim = size(Im);

p_n1 = zeros(dim(1), dim(2), 2, 3);

u_n1 = zeros(dim(1), dim(2), 3);

dif = inf;
Iter = 0;
while dif > r && Iter < maxIter
    
    
    Iter = Iter + 1;
    u_n = u_n1;
    f = Im.*(1 - mask) + u_n1.*mask;
    p_n = p_n1;
    norm_p = sum(p_n(:).^2);
    
    if  norm_p > 1
       p_n = p_n/norm_p;
    end    
    
    %Compute p_n+1
    div_p_1 = Divergence_p( p_n(: ,: , 1, 1 ), p_n(: ,: , 2, 1) );
    [ grad_v_x, grad_v_y ] = Gradient_v( div_p_1 - (f( :, :, 1)/lambda) );
    grd_div_1 = zeros(dim(1), dim(2), 2);
    grd_div_1(:, :, 1) =  grad_v_x;
    grd_div_1(:, :, 2) =  grad_v_y ;
    
    div_p_2 = Divergence_p( p_n(: ,: , 1, 2 ), p_n(: ,: , 2, 2) );
    [ grad_v_x, grad_v_y ] = Gradient_v( div_p_2 - (f(:, :, 2)/lambda) );
    grd_div_2 = zeros(dim(1), dim(2), 2);
    grd_div_2(:, :, 1) =  grad_v_x;
    grd_div_2(:, :, 2) =  grad_v_y ;
    
    div_p_3 = Divergence_p( p_n(: ,: , 1, 3 ), p_n(: ,: , 2, 3) );
    [ grad_v_x, grad_v_y ] = Gradient_v( div_p_3 - (f(:, :, 3)/lambda) );
    grd_div_3 = zeros(dim(1), dim(2), 2);
    grd_div_3(:, :, 1) =  grad_v_x;
    grd_div_3(:, :, 2) =  grad_v_y ;
    
    % i = 1
    p_n1(:, :, :, 1) = (p_n(:, :, :, 1 ) + delta_t*grd_div_1)./...
        ( 1 + delta_t*sqrt(grd_div_1.^2 + grd_div_2.^2 + grd_div_3.^2));
    % i = 2
    p_n1(:, :, :, 2) = (p_n(:, :, :, 2 ) + delta_t*grd_div_2)./...
        ( 1 + delta_t*sqrt(grd_div_1.^2 + grd_div_2.^2 + grd_div_3.^2));
    % i = 3
    p_n1(:, :, :, 3) = (p_n(:, :, :, 3 ) + delta_t*grd_div_3)./...
        ( 1 + delta_t*sqrt(grd_div_1.^2 + grd_div_2.^2 + grd_div_3.^2));
    
    div_p_1 = Divergence_p( p_n1(: ,: , 1, 1 ), p_n1(: ,: , 2, 1) );
    div_p_2 = Divergence_p( p_n1(: ,: , 1, 2 ), p_n1(: ,: , 2, 2) );
    div_p_3 = Divergence_p( p_n1(: ,: , 1, 3 ), p_n1(: ,: , 2, 3) );
    
    u_n1(:, :, 1) = f(:, :, 1) - lambda*div_p_1;
    u_n1(:, :, 2) = f(:, :, 2) - lambda*div_p_2;
    u_n1(:, :, 3) = f(:, :, 3) - lambda*div_p_3;
 
    dif = max(max(max(abs(u_n1 - u_n))));
    
end 
f = Im.*(1 - mask) + u_n1.*mask;
Iinp = f;
end

