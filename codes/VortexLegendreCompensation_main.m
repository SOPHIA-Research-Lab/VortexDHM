%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:-->        VortexLegendreCompensation_main
%
%
% Description:-->  This algorithm allows to evaluate the performance of phase
%                  compensation from two steps: 1) Solving the Vortex subpixel
%                  frequencies and 2) fitting square Legendre polynomials.
%
%
% Authors:-->      Rene Restrepo(2), Karina Ortega-Sánchez(1,2), Carlos Trujillo(2)
%                  Ana Doblas(3) Alfonso Padilla (1)
%                  (1) Universidad Politecnica de Tulancingo
%                  (2) EAFIT Univeristy SOPHIA Research Group (Applied Optics)
%                  (3) Department of Electrical and Computer Engineering, University 
%                      of Massachusetts Dartmouth
% Email:-->        rrestre6@eafit.edu.co
% Date:-->         12/17/2024
% version 1.0 (2024)
% Notes-->


%% clear memory, worksapce,close all figures and clear command window (memory variables)
clc % clean
close all% close all windows
clear all% clear of memory all variable

%% Lines to add folders for reading images and/or functions
addpath('../vortex/HoloPaper');

%% Load the hologram

name =  'NredBlood_40x_632_4.65_30mm.tiff';
% name =  'TDrosophila_0x_632_5.33.bmp';
% name =  'Tesper_0x_632_3.75.bmp';
% name =  'Tglibocut_40x_532_3.45.jpg';
% name =  'TProbiotics_20x_632_3.75.bmp';
% name =  'Tslm1_10x40x_632_3.75.bmp';
% name =  'NStar_20x_532_5.86_-4cm.tiff';

hologram = double(imread(name));
[~, nameNoExt] = fileparts(name);
%% Verified that the hologram has a square size
hologram = hologram(:,:,1);
if size(hologram,2) ~= size(hologram,1)
    lim=round(size(hologram,2)-size(hologram,1))/2;
    hologram = hologram(:,lim:end-lim-1);
end

[N, M]=size(hologram);
[n,m] = meshgrid(-M/2:M/2-1,-N/2:N/2-1);
%% Lines to load parameters and define constants
posUnderscore = strfind(nameNoExt, '_');

if length(posUnderscore) < 3
    lambda = 0.632;
    dxy = 3.75;
else
    lambda_str = (nameNoExt(posUnderscore(2) + 1 : posUnderscore(2) + 3));
    lambda = str2double(['0.' lambda_str]);
    dxy = str2double(nameNoExt(posUnderscore(3) + 1 : posUnderscore(3) + 4));
end
k = 2 * pi / lambda;
fx_0 = M/2;
fy_0 = N/2;

%% Algorithm parameters to change
medFilt = 3; % Change between 3 or 1
NoPistonCompensation = false; % "False" to enable piston compensation.
UsePCA = true; % "True" to enable compensation using PCA.
limit= 64/2; % Spatial frequency support region region for calculating the Legendre coefficients.

%% Filter oreder +1
[holo_filtered, fxOverMax, fyOverMax, cir_mask] = ...
    functions_evaluation.spatial_filter(hologram,M,N,'0',3);
%% Vortex subpixel frequencies calculation

ft_holo = log(abs(fftshift(fft2(fftshift(hologram)))).^2);
field = medfilt2(ft_holo, [medFilt, medFilt], 'symmetric');
[positions] = vortexCompensation(field, fxOverMax, fyOverMax);


%% Tilt compensation aberration using optical vortex

[ref_wave] = functions_evaluation.reference_wave...
    (M,N,m,n,lambda,dxy,positions(1),positions(2),k,fx_0,fy_0);
field_compensate = ref_wave.*holo_filtered;

figure,imagesc(angle(field_compensate)),colormap(gray),colorbar
title('Phase compensated (Vortex)'),daspect([1 1 1])

[phaseCorrected, Legendre_Coefficients] = ...
    LegendreCompensation(field_compensate, limit, NoPistonCompensation, UsePCA);

%% Vortex + piston

if ~NoPistonCompensation
    PistonCorrtection = ones(size(field_compensate)).*Legendre_Coefficients(1);
    vortexPistonCorrection = field_compensate.*(exp(1i.*PistonCorrtection));
    figure,imagesc(angle(vortexPistonCorrection)),colormap(gray),colorbar
    title('Phase compensated (Vortex + Piston)'),daspect([1 1 1])
end
gridSize = size(angle(field_compensate),1);
[X, Y] =  meshgrid(-1:(2 / gridSize):(1 - 2 / gridSize), ...
    -1:(2 / gridSize):(1 - 2 / gridSize));
dA=(2 / gridSize) ^ 2;

order = 2:6;
[polynomials] = squareLegendrefitting(order, X, Y);
Legendres = reshape(polynomials, [size(polynomials, 1)*size(polynomials, 2) ...
    size(polynomials, 3)]);
zProds = Legendres.'* Legendres * dA;
zNorm = bsxfun(@rdivide, Legendres, sqrt(diag(zProds).'));
Legendres = (bsxfun(@times, (ones(size(order))').', zNorm));
Legendres_norm_const =sum(Legendres.^2,1)*dA;

WavefrontReconstructed_Vect = sum(repmat(Legendre_Coefficients(2:size(order,2)+1)./ ...
    sqrt(Legendres_norm_const(1:size(order,2))),[size(Legendres,1) 1]).*Legendres(:,1:size(order,2)),2);

WavefrontReconstructed = reshape(WavefrontReconstructed_Vect, ...
    [size(polynomials, 1) size(polynomials, 2)]);

compensatedHologram = (abs(field_compensate).*exp(1i.* angle(field_compensate))./ ...
            exp((1i).*WavefrontReconstructed));

figure,imagesc(angle(compensatedHologram)),colormap(gray),colorbar
title('Phase compensated (Vortex + Legendre)'),daspect([1 1 1])

if ~NoPistonCompensation
    PistonCorrtection = ones(size(compensatedHologram)).*Legendre_Coefficients(1);
    vortexPistonCorrection = compensatedHologram.*(exp(1i.*PistonCorrtection));
    figure,imagesc((angle(vortexPistonCorrection))),colormap(gray),colorbar
    title('Phase compensated (Vortex + Legendre + Piston)'),daspect([1 1 1])

end
