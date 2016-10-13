clear all
I=imread('image1_toRestore.jpg');
pic=double(I);
[m n]= size(I);
[x,y] = meshgrid(1:n, 1:m); % rejilla imagen inicial
r=0.1; % factor de escala
[p,q]=meshgrid(1:r:n, 1:r:m); % rejilla imagen final
I2=interp2(x,y,I,p,q,'bicubic'); % interpolación
% 'nearest','bilinear','bicubic'
figure
subplot(1,2,1),imagesc(I),axis image
title('Original','FontSize',18)
subplot(1,2,2),imagesc(I2),axis image
title('Interpolador NN ','FontSize',18)
