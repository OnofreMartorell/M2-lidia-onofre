%close all;
clearvars;
clc
addpath(genpath('../Block2'))

% I = double(imread('circles.png'));
% I = double(imread('noisedCircles.tif')); %Lambda 7, 6.5
% I = double(imread('phantom17.bmp')); %Lambda 6.5
% I = double(imread('phantom18.bmp')); %Lambda 12
I = double(imread('phantom19.bmp'));


% I = double(imread('Tres_formas.png'));

I = mean(I,3);
I = I - min(I(:));
I = I/max(I(:));

[ni, nj] = size(I);

params_ROF.r = 10^-4;
params_ROF.maxIter = 1000;
params_ROF.delta_t = 1/4;


%%Parameters
param.lambda = 12;
%lambda1=10^-3; %Hola carola problem
%lambda2=10^-3; %Hola carola problem


param.tol = 1000;
param.iterMax = 1;

[X, Y] = meshgrid(1:nj, 1:ni);

%%Initial phi
% E_0 = (-sqrt( ( X-round(ni/2)).^2 + (Y-round(nj/2)).^2) + 50);
%%% This initialization allows a faster convergence for phantom 18
E_0 = (-sqrt( ( X-round(ni/2)).^2 + (Y-round(nj/4)).^2)+50);

E_0 = E_0 > 0;

%phi_0=I; %For the Hola carola problem
% Normalization of the initial phi to [-1 1]
% E_0 = E_0 - min(E_0(:));
% E_0 = 2*E_0/max(E_0(:));
% E_0 = E_0 - 1;


%%Explicit Gradient Descent
seg = segmentation_Chambolle( I, E_0, param, params_ROF );
sprintf('Finished')