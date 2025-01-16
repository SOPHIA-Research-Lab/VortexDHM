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
%                  Ana Doblas(2) Alfonso Padilla (1)
%                  (1) Universidad Politecnica de Tulancingo
%                  (2) EAFIT Univeristy SOPHIA Research Group (Applied Optics)                                                                                                   
% Email:-->        rrestre6@eafit.edu.co                                         
% Date:-->         12/17/2024                                                      
% version 1.0 (2024)                                                               
% Notes-->                                                                         


%% clear memory, worksapce,close all figures and clear command window (memory variables)
clc % clean
close all% close all windows
clear all% clear of memory all variable

%% Lines to add folders for reading images and/or functions
addpath('../holograms');

%% Load the hologram 
name = 'starTarget.bmp'; %StarTraget telecentric
% name =  'usaf.bmp'; %USAF telecentric
% name ='-4cm_20x_star.tiff'; %Star target no telecentric
% name =  'redBlood30mm.tiff'; % Red blood cell no telecentric

[hologram,M,N,m,n] = functions_evaluation.holo_read(name);

%% Lines to load parameters and define constants
%Use these parameters for the starTarget.bmp and the usaf.bmp
lambda = 0.632;
dxy = 3.75;

% Use these parameters for the redBlood30mm.tiff and the -4cm_20x_star.tiff
% lambda = 0.532;
% dxy = 4.75;
k = 2 * pi / lambda;
fx_0 = M/2;
fy_0 = N/2;

%% Vortex subpixel frequencies calculation
% Selection order (+1 DO)
[holo_filtered, fxOverMax, fyOverMax, cir_mask] = ...
    functions_evaluation.spatial_filter(hologram,M,N,'Yes',3);

% Tilt compensation aberration using by optical vortex
[positions] = vortexCompensation(hologram, cir_mask, fxOverMax, fyOverMax);

%% Hologrom correction
[ref_wave] = functions_evaluation.reference_wave...
    (M,N,m,n,lambda,dxy,positions(1),positions(2),k,fx_0,fy_0);

field_compensate = ref_wave.*holo_filtered;
figure,imagesc(angle(field_compensate)),colormap(gray),colorbar
title('Phase compensated'),daspect([1 1 1])

%% Square Legendre polynomials fitting
gridSize = size(angle(field_compensate),1);
if size(field_compensate,2) ~= size(field_compensate,1)
    limit=round(size(field_compensate,2)-size(field_compensate,1))/2;
    field_compensate = field_compensate(:,limit:end-limit-1);
end

[X, Y] =  meshgrid(-1:(2 / gridSize):(1 - 2 / gridSize), ...
    -1:(2 / gridSize):(1 - 2 / gridSize));
dA=(2 / gridSize) ^ 2;
order = 1:5;
[polynomials] = squareLegendrefitting(order, X, Y);

% Orthonormalization of Legendre polynomials
Legendres = reshape(polynomials, [size(polynomials, 1)*size(polynomials, 2) ...
    size(polynomials, 3)]);
zProds = Legendres.'* Legendres * dA;
zNorm = bsxfun(@rdivide, Legendres, sqrt(diag(zProds).'));
Legendres = (bsxfun(@times, (ones(size(order))').', zNorm));
Legendres_norm_const =sum(Legendres.^2,1)*dA;

phaseVector = reshape(phase_unwrap(angle(field_compensate(1:gridSize, 1:gridSize))), ...
    [size(polynomials, 1)*size(polynomials, 2) 1]);

% Coefficients calculation of Legendre polynomials
Legendre_Coefficients=sum(bsxfun(@times, Legendres, phaseVector), 1)*dA;

WavefrontReconstructed_Vect = sum(repmat(Legendre_Coefficients(1:end)./ ...
    sqrt(Legendres_norm_const(1:end)),[size(Legendres,1) 1]).*Legendres(:,1:end),2);
WavefrontReconstructed = reshape(WavefrontReconstructed_Vect, [size(polynomials, 1) size(polynomials, 2)]);

PhaseNoPistonAndTilts = (exp(1i.* angle(field_compensate(1:gridSize, 1:gridSize)))./ ...
    exp((1i).*WavefrontReconstructed));

figure(4), imagesc(angle(PhaseNoPistonAndTilts)), colorbar, colormap gray


