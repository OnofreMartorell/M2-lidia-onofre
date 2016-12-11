function [ div ] = Divergence_p( p_i_x, p_i_y )
[N_x, N_y] = size(p_i_x);
div = zeros(N_x, N_y);

for i_x = 1:N_x
    for i_y = 1:N_y
        if i_x == 1
            if i_y == 1
                div(i_x, i_y) = p_i_x(i_x, i_y) + p_i_y(i_x, i_y);
            else
                if i_y == N_y
                    div(i_x, i_y) = p_i_x(i_x, i_y) - p_i_y(i_x, i_y - 1);
                else
                    div(i_x, i_y) = p_i_x(i_x, i_y) + p_i_y(i_x, i_y) - p_i_y(i_x, i_y - 1);
                end
            end
        else
            if i_x == N_x
                if i_y == 1
                    div(i_x, i_y) = - p_i_x(i_x - 1, i_y) + p_i_y(i_x, i_y);
                else
                    if i_y == N_y
                        div(i_x, i_y) = - p_i_x(i_x - 1, i_y) - p_i_y(i_x, i_y - 1);
                    else
                        div(i_x, i_y) = - p_i_x(i_x - 1, i_y) + p_i_y(i_x, i_y) - p_i_y(i_x, i_y - 1);
                    end
                end
            else
                if i_y == 1
                    div(i_x, i_y) = p_i_x(i_x, i_y) - p_i_x(i_x - 1, i_y) + p_i_y(i_x, i_y);
                else
                    if i_y == N_y
                        div(i_x, i_y) = p_i_x(i_x, i_y) - p_i_x(i_x - 1, i_y) - p_i_y(i_x, i_y - 1);
                    else
                        div(i_x, i_y) = p_i_x(i_x, i_y) - p_i_x(i_x - 1, i_y) + p_i_y(i_x, i_y) - p_i_y(i_x, i_y - 1);
                    end
                end
            end
        end
    end
end

