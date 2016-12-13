function [ Iinp ] = G7_inpainting_color( Im, mask, param )
lambda = param.lambda;
delta_t = param.delta_t;

r = param.r;

maxIter = param.maxIter;


dim = size(Im);

p_n1 = zeros(dim(1), dim(2), 2, 3);

u_n1 = zeros(dim(1), dim(2), 3);

dif = inf;
Iter = 0;
while dif > r && Iter < maxIter
    
    
    Iter = Iter + 1;
    Iter
    
    u_n = u_n1;
    f = Im.*(1 - mask) + u_n1.*mask;
    p_n = p_n1;
    
    %Compute p_n+1
    grd_div = zeros(dim(1), dim(2), 2, 3);
    div_p = zeros(dim(1), dim(2), 3);
    for i = 1:3
        div_p(:, :, i) = Divergence_p( p_n(: ,: , 1, i ), p_n(: ,: , 2, i) );
        [ grad_v_x, grad_v_y ] = Gradient_v( div_p(:, :, i) - (f( :, :, i)/lambda) );
        grd_div(:, :, 1, i) =  grad_v_x;
        grd_div(:, :, 2, i) =  grad_v_y ;
    end
    
    for i = 1:3
        p_n1(:, :, :, i) = (p_n(:, :, :, i ) + delta_t*grd_div(:, :, :, i))./...
            ( 1 + delta_t*sqrt(grd_div(:, :, :, 1).^2 + grd_div(:, :, :, 2).^2 + grd_div(:, :, :, 3).^2));
        
    end
    norm_p1 = sqrt(sum(p_n1(:).^2));
    
    if  norm_p1 > 1
        p_n1 = p_n1/norm_p1;
    end
    
    for i = 1:3
        div_p(:, :, i) = Divergence_p( p_n1(: ,: , 1, i ), p_n1(: ,: , 2, i) );
        [ grad_v_x, grad_v_y ] = Gradient_v( div_p(:, :, i) - (f( :, :, i)/lambda) );
        grd_div(:, :, 1, i) =  grad_v_x;
        grd_div(:, :, 2, i) =  grad_v_y ;
    end
    
    for i = 1:3
        u_n1(:, :, i) = f(:, :, i) - lambda*div_p(:, :, i);
    end
    
    
    dif = max(max(max(abs(u_n1 - u_n))));
%     figure(3)
%     imshow(f)
end

f = Im.*(1 - mask) + u_n1.*mask;
Iinp = f;
end
