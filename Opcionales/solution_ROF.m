function [ u ] = solution_ROF( Im, param )
lambda = param.lambda;
delta_t = param.delta_t;

r = param.r;

maxIter = param.maxIter;

[ni, nj, nC] = size(Im);
dim = [ni nj];
f = Im;
p_n1 = zeros(dim(1), dim(2), 2, nC);

u_n1 = zeros(dim(1), dim(2), nC);

dif = inf;
Iter = 0;
while dif > r && Iter < maxIter
    
    Iter = Iter + 1;
    
    u_n = u_n1;
    p_n = p_n1;
    
    %Compute p_n+1
    grd_div = zeros(dim(1), dim(2), 2, nC);
    div_p = zeros(dim(1), dim(2), nC);
    
    for i = 1:nC
        div_p(:, :, i) = Divergence_p( p_n(: ,: , 1, i ), p_n(: ,: , 2, i) );
        [ grad_v_x, grad_v_y ] = Gradient_v( div_p(:, :, i) - (f( :, :, i)/lambda) );
        grd_div(:, :, 1, i) =  grad_v_x;
        grd_div(:, :, 2, i) =  grad_v_y ;
    end
    
    for i = 1:nC
        p_n1( :, :, :, i) = (p_n( :, :, :, i ) + delta_t*grd_div(:, :, :, i))./...
            ( 1 + delta_t*sqrt(sum(grd_div.^2, 4)));
        
    end
    
    norm_p1 = sqrt(sum(p_n1(:).^2));
    
    if  norm_p1 > 1
        p_n1 = p_n1/norm_p1;
    end
    
    
    for i = 1:nC
        u_n1(:, :, i) = f(:, :, i) - lambda*div_p(:, :, i);
    end
    
    
    dif = max(max(max(abs(u_n1 - u_n))));
%     figure(3)
%     imshow(f)
end

u = u_n1;

end

