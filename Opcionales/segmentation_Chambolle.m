function [ seg ] = segmentation_Chambolle( I, E_0, params_seg, params_ROF )



E = E_0;
dif = inf;
nIter = 0;

figure(1)
imagesc(I);
colormap gray;
hold on;
E_compl = E == 0;
partition = double(E) - double(E_compl);
contour(partition, [0 0], 'r', 'LineWidth', 2)
title('Segmentation over the image')

axis off;
hold off
drawnow;

while dif > params_seg.tol && nIter < params_seg.iterMax
    
    E_old = E;
    nIter = nIter + 1;
    nIter
    
    %Minimization w.r.t c1 and c2 (constant estimation)
    c1 = sum(sum(I.*E_old))/(sum(sum(E_old))); 
    c2 = sum(sum(I.*(1 - E_old)))/(sum(sum(1 - E_old))); 

    params_ROF.lambda = params_seg.lambda/(c2 - c1);
    u = solution_ROF( I, params_ROF );
    E = u > (c1 + c2)/2;
    
    dif = sum( (E(:) - E_old(:)).^2 );
    
    %Plot the curve evolution over the image
    figure(1)
    imagesc(I);
    colormap gray;
    hold on;
    E_compl = E == 0;
    partition = double(E) - double(E_compl);
    contour(partition, [0 0], 'r', 'LineWidth', 2) 
    title('Segmentation over the image')
    
    axis off;
    hold off
    drawnow;
    
end
seg = E;
end

