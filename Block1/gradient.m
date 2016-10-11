function [ grad_v_x, grad_v_y ] = gradient( v )
[N_x, N_y]  = size(v);
grad_v_x = zeros(N_x, N_y);
grad_v_y = zeros(N_x, N_y);

for i_x = 1:N_x
    for i_y = 1:N_y
        if i_x  == N_x
            grad_v_x(i_x, i_y) = 0;
        else
            grad_v_x(i_x, i_y) = v(i_x + 1, i_y) - v(i_x, i_y);
        end  
        
        if i_y  == N_y
            grad_v_y(i_x, i_y) = 0;
        else
            grad_v_y(i_x, i_y) = v(i_x, i_y + 1) - v(i_x, i_y);
        end
    end
end
end

