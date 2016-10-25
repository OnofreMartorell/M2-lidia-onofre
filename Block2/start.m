%close all;
clearvars;
clc

%I=double(imread('zigzag_mask.png'));
%I=mean(I,3); %To 2D matrix
% I = double(imread('circles.png'));
%I=double(imread('noisedCircles.tif'));
% I=double(imread('phantom17.bmp'));
%I=double(imread('phantom18.bmp'));

I = double(imread('Tres_formas.png'));

I = mean(I,3);
I = I - min(I(:));
I = I/max(I(:));

[ni, nj] = size(I);



%Lenght and area parameters
%circles.png mu=1, mu=2, mu=10
%noisedCircles.tif mu=0.1
%phantom17 mu=1, mu=2, mu=10
%phantom18 mu=0.2 mu=0.5
%hola carola mu=1
mu = 1;
nu = 0;


%%Parameters
lambda1 = 1;
lambda2 = 1;
%lambda1=10^-3; %Hola carola problem
%lambda2=10^-3; %Hola carola problem

epHeaviside = 1;
%eta=0.01;

eta = 1;
tol = 0.00000000001;
%dt=(10^-2)/mu; 
dt = (10^-1)/mu;
iterMax = 1000;
%reIni=0; %Try both of them
%reIni=500;

reIni = 1000;
[X, Y] = meshgrid(1:nj, 1:ni);

%%Initial phi
phi_0 = (-sqrt( ( X-round(ni/2)).^2 + (Y-round(nj/2)).^2)+50);
% phi_0 = sin(pi/5*X) * sin(pi/5*Y);
%%% This initialization allows a faster convergence for phantom 18
%phi_0=(-sqrt( ( X-round(ni/2)).^2 + (Y-round(nj/4)).^2)+50);



%phi_0=I; %For the Hola carola problem
% Normalization of the initial phi to [-1 1]
phi_0 = phi_0 - min(phi_0(:));
phi_0 = 2*phi_0/max(phi_0(:));
phi_0 = phi_0 - 1;




%%Explicit Gradient Descent
seg = G7_sol_ChanVeseIpol_GDExp( I, phi_0, mu, nu, eta, lambda1, lambda2, tol, epHeaviside, dt, iterMax, reIni );

%% Goal Image
% Computation
clearvars;

%Read the image
Im = double(imread('image_to_Restore.png'));

[ni, nj, ~] = size(Im);
Im = Im - min(Im(:));
Im = Im/max(Im(:));


%We want to inpaint those areas in which mask == 1 (red part of the image)
I_ch1 = Im(:,:,1);
I_ch2 = Im(:,:,2);
I_ch3 = Im(:,:,3);


%Create the mask
mask = zeros(ni,nj) ; 
%mask_img(i,j) == 1 means we have lost information in that pixel
                                      %mask(i,j) == 0 means we have information in that pixel
for i = 1:ni
    for j = 1:nj
        if(I_ch1(i,j) == 1 && I_ch2(i,j) == 0 && I_ch3(i,j) == 0 && i > 284)
            mask(i,j) = 1;
        end
    end
end


I = mask;

%Lenght and area parameters
%circles.png mu=1, mu=2, mu=10
%noisedCircles.tif mu=0.1
%phantom17 mu=1, mu=2, mu=10
%phantom18 mu=0.2 mu=0.5
%hola carola mu=1
mu = 1;
nu = 0;


%%Parameters
% lambda1 = 1;
% lambda2 = 1;
lambda1 = 10^-2; %Hola carola problem
lambda2 = 10^-4; %Hola carola problem

epHeaviside = 1;
%eta=0.01;

eta = 1;
tol = 0.0001;
%dt=(10^-2)/mu; 
dt = (10^-1)/mu;
iterMax = 2500;
%reIni=0; %Try both of them
%reIni=500;

reIni = 1000;
[X, Y] = meshgrid(1:nj, 1:ni);

%%Initial phi
% phi_0 = (-sqrt( ( X-round(ni/2)).^2 + (Y-round(nj/2)).^2)+50);
% phi_0 = sin(pi/5*X) * sin(pi/5*Y);
%%% This initialization allows a faster convergence for phantom 18
%phi_0=(-sqrt( ( X-round(ni/2)).^2 + (Y-round(nj/4)).^2)+50);
%For the Hola carola problem
phi_0 = I; 


% Normalization of the initial phi to [-1 1]
phi_0 = phi_0 - min(phi_0(:));
phi_0 = 2*phi_0/max(phi_0(:));
phi_0 = phi_0 - 1;

seg2 = G7_sol_ChanVeseIpol_GDExp( I, phi_0, mu, nu, eta, lambda1, lambda2, tol, epHeaviside, dt, iterMax, reIni );