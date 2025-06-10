%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:-->        Legendre compensation over a reduced square area                       
%                                                                                                               
%                                                                                  
% Description:-->  This algorithm calculates the principal components over a reduced 
%                  region of size m x n and computes the Legendre coefficients.           
%                                                                                  
% Authors:-->      Rene Restrepo(2), Karina Ortega-Sanchez(1,2), 
%                  Carlos Trujillo(2), Ana Doblas(3) and Alfonso Padilla (1)
%                  (1) Polytechnic University of Tulancingo
%                  (2) EAFIT Univeristy SOPHIA Research Group (Applied Optics) 
%                  (3) Department of Electrical and Computer Engineering, University 
%                      of Massachusetts Dartmouth
%
% Inputs:-->       (1) Field with tilt compensation 
%                  (2) Limit that defines the size of the reduced region
%                  (3) Piston compensation that can take "true" or "false" values
%                  (4) Principal component calculation that can take "true" or "false" values
% Outputs:-->      (1) Phase corrected
%                  (2) Legendre Coefficients
%
% Email:-->        rrestre6@eafit.edu.co                                         
% Date:-->         12/17/2024                                                      
% version 1.0 (2024)                                                               
% Notes-->    

function [compensatedHologram, Legendre_Coefficients] ...
    = LegendreCompensation(field_compensate, limit, NoPistonCompensation, UsePCA)
fftFieldCompensated = fftshift(fft2(ifftshift(field_compensate)));

%% Reduced region calculation
[A, B] = size(fftFieldCompensated);
fftFieldCompensated = ...
    fftFieldCompensated(round(A/2)+1-limit:round(A/2)+limit, round(B/2)+1-limit:round(B/2)+limit);
square = ifftshift(ifft2(fftshift(fftFieldCompensated)));

%% Calculation of the first dominant component using Singular Value Decomposition (SVD)
if UsePCA
    [U, S, V] = svd((square));
    numComponentes = 1;
    dominantComponent = U(:, 1:numComponentes) * ...
        S(1:numComponentes, 1:numComponentes) * V(:, 1:numComponentes)';
    dominantComponent = phase_unwrap(angle(dominantComponent));

else
    dominantComponent = phase_unwrap(angle(square));
end

%% Specify the number of Legendre coefficients to calculate
gridSize = size(dominantComponent,1);
[X, Y] =  meshgrid(-1:(2 / gridSize):(1 - 2 / gridSize), ...
    -1:(2 / gridSize):(1 - 2 / gridSize));
dA=(2 / gridSize) ^ 2;
order = 1:10;
[polynomials] = squareLegendrefitting(order, X, Y);

%% Orthonormalization of Legendre polynomials
Legendres = reshape(polynomials, [size(polynomials, 1)*size(polynomials, 2) ...
    size(polynomials, 3)]);
zProds = Legendres.'* Legendres * dA;
zNorm = bsxfun(@rdivide, Legendres, sqrt(diag(zProds).'));
Legendres = (bsxfun(@times, (ones(size(order))').', zNorm));
Legendres_norm_const =sum(Legendres.^2,1)*dA;

phaseVector = reshape(dominantComponent, ...
    [size(polynomials, 1)*size(polynomials, 2) 1]);
%% Calculation of Legendre coefficients
Legendre_Coefficients=sum(bsxfun(@times, Legendres, phaseVector), 1)*dA;

%% Piston correction 
if NoPistonCompensation
    WavefrontReconstructed_Vect = sum(repmat(Legendre_Coefficients(2:end)./ ...
        sqrt(Legendres_norm_const(2:end)),[size(Legendres,1) 1]).*Legendres(:,2:end),2);

    WavefrontReconstructed = reshape(WavefrontReconstructed_Vect, ...
        [size(polynomials, 1) size(polynomials, 2)]);

    compensatedHologram = (exp(1i.* angle(square))./ ...
        exp((1i).*WavefrontReconstructed));
else
    evaluationPiston = 1;
    sizePistonEvaluation = -pi:(pi/6):pi;
    for i = sizePistonEvaluation
        Legendre_Coefficients(1) = i;
        WavefrontReconstructed_Vect = sum(repmat(Legendre_Coefficients(1:end)./ ...
            sqrt(Legendres_norm_const(1:end)),[size(Legendres,1) 1]).*Legendres(:,1:end),2);
        WavefrontReconstructed = reshape(WavefrontReconstructed_Vect, ...
            [size(polynomials, 1) size(polynomials, 2)]);
        compensatedHologram = (exp(1i.* angle(square))./ ...
            exp((1i).*WavefrontReconstructed));

        evaluationPistonPhase = angle(compensatedHologram);
        variance(evaluationPiston)=var(evaluationPistonPhase(:));
        evaluationPiston = evaluationPiston+1;
    end
%% Determining the piston value for proper compensation using variance analysis
    [minVz, posVz] = min(variance);
    Legendre_Coefficients (1) = sizePistonEvaluation(posVz);
    WavefrontReconstructed_Vect = sum(repmat(Legendre_Coefficients(1:end)./ ...
        sqrt(Legendres_norm_const(1:end)),[size(Legendres,1) 1]).*Legendres(:,1:end),2);
    WavefrontReconstructed = reshape(WavefrontReconstructed_Vect, ...
        [size(polynomials, 1) size(polynomials, 2)]);
    compensatedHologram = (exp(1i.* angle(square))./ ...
        exp((1i).*WavefrontReconstructed));
end






