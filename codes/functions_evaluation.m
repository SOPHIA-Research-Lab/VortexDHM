%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: functions_evaluation                                                   %
%                                                                              %
% The script contains all implemented function for the evaluation_main.pyget   %
%                                                                              %                                       
% Authors: Raul Castaneda and Ana Doblas                                       %
% Applied Optics Group EAFIT univeristy                                        % 
%                                                                              %
% Email: racastaneq@eafit.edu.co; adoblas@umassd.edu                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef functions_evaluation
    methods(Static)

        function [holo,M,N,m,n] = holo_read(filename)
            holo = double(imread(filename));
            holo = holo(:,:,1);
            [N,M] = size(holo);
            [n,m] = meshgrid(-M/2:M/2-1,-N/2:N/2-1);
        end


        function [holo_filtered,fx_max,fy_max,cir_mask] = spatial_filter(holo,M,N,visual,factor)
            % Compute Fourier Transfor hologram
            ft_holo = fftshift(fft2(fftshift(holo)));
            ft_holo(1:5,1:5)=0;
            % filter DC term
            mask1 = ones(N,M);
            mask1(N/2-20:N/2+20,M/2-20:M/2+20)=0;
            ft_holo_I = ft_holo .* mask1;
        
            % filter reflection
            mask1 = ones(N,M);
            mask1(1,1)=0;
            ft_holo_I = ft_holo_I .* mask1;
            
            % select region of interest (ROI) to search the diffraction term.
            region_interest = ft_holo_I(1:N, 1:M/2);
            [maxValue, linearIndex] = max(abs(region_interest), [], 'all', 'linear');
            [fy_max, fx_max] = ind2sub(size(region_interest), linearIndex);
            
            % ciruclar mask to filter the +1 difraction term
            distance = sqrt((fx_max - N/2)^2+(fy_max - M/2)^2);
            resc = distance/factor;
            cir_mask = ones(N,M);
            for r=1:N
                for p=1:M
                    if sqrt((r-fy_max)^2+(p-fx_max)^2)>resc
                        cir_mask(r,p)=0;
                    end
                end
            end
    
            % Applying spatical filter and computing the filtered hologram
            ft_holo_filtered = ft_holo .* cir_mask;
            holo_filtered = fftshift(ifft2(fftshift(ft_holo_filtered)));
            
            if visual == 'Yes'
                figure,imagesc(log(abs(ft_holo).^2)),colormap(gray),title('FT Hologram'),daspect([1 1 1])
                figure,imagesc(cir_mask),colormap(gray),title('Circular Filter'),daspect([1 1 1])
                figure,imagesc(log(abs(ft_holo_filtered).^2)),colormap(gray),title('FT Filtered Hologram'),daspect([1 1 1])
                saveas(gcf, 'mask.bmp');
            end
        end


        function [ref_wave] = reference_wave(M,N,m,n,lambda,dxy,fx_max,fy_max,k,fx_0,fy_0)
            theta_x = asin((fx_0 - fx_max) * lambda / (M * dxy));
            theta_y = asin((fy_0 - fy_max) * lambda / (N * dxy));
            ref_wave = exp(1i * k * (sin(theta_x) * n * dxy + sin(theta_y) * m * dxy));
        end

        function [cf] = minimization_LegCoeff(seed, gridSize, Legendres, Legendre_Coefficients, Legendres_norm_const, polynomials, field_compensate)
            cf = 0;

            Legendre_Coefficients=[seed Legendre_Coefficients(1:end)];
            Legendres_norm_const=[1 Legendres_norm_const(1:end)];
            Legendres=[ones(length(Legendres),1) Legendres];

            WavefrontReconstructed_Vect = sum(repmat(Legendre_Coefficients(1:end)./ ...
                                        sqrt(Legendres_norm_const(1:end)),[size(Legendres,1) 1]).*Legendres(:,1:end),2);
            WavefrontReconstructed = reshape(WavefrontReconstructed_Vect, [size(polynomials, 1) size(polynomials, 2)]);

            PhaseNoPistonAndTilts = (exp(1i.* angle(field_compensate(1:gridSize, 1:gridSize)))./ ...
                                    exp((1i).*WavefrontReconstructed));

            phase = angle(PhaseNoPistonAndTilts);
            [M,N] = size(phase);
            ib = imbinarize(phase, 0.5);
            cf = M*N - sum(ib(:));
            
        end


        function [holo_filter,holo_FT,fx_max,fy_max] = spatialFilter_SIDHM(holo,M,N,visual,factor)
            fft_holo = fftshift(fft2(fftshift(holo(:,:,1)))); 
            figure,imagesc(log(abs(fft_holo).^2)),colormap(gray),title('FT Hologram'),daspect([1 1 1]) 
            
            [pointsY,pointsX] = ginput(2);
            x3 = round(pointsX(1));
            x4 = round(pointsX(2));
            y3 = round(pointsY(1));
            y4 = round(pointsY(2));
           
            mask = zeros(M,N);
            mask(x3:x4,y3:y4) = 1;
            fft_holo_I = fft_holo .* mask;
            figure,imagesc((abs( fft_holo_I).^0.1)),colormap(gray),title('FT Hologram'),daspect([1 1 1]) 
            
            % max values first peak
            maxValue_1 = max(max(abs(fft_holo_I)));
            [fy_max_1, fx_max_1] = find(abs(fft_holo_I) == maxValue_1);
            mask(fy_max_1 - 20:fy_max_1 + 20,fx_max_1 - 20:fx_max_1 + 20)=0;
            fx_max_L = fx_max_1;
            fy_max_L = fy_max_1;
            fft_holo_I = fft_holo_I .* mask;
            figure,imagesc(log(abs(fft_holo_I).^2)),colormap(gray),title('FT FIrst peak'),daspect([1 1 1])
            
            maxValue_1 = max(max(abs(fft_holo_I)));
            [fy_max_1, fx_max_1] = find(abs(fft_holo_I) == maxValue_1);
            fx_max_D = fx_max_1(1);
            fy_max_D = fy_max_1(1);
        
            fx_max = [fx_max_L,fx_max_D];
            fy_max = [fy_max_L,fy_max_D];
        
            %find the centers between both peaks 
            middlePoint_X = (fx_max_D + fx_max_L) / 2;
            middlePoint_Y = (fy_max_D + fy_max_L) / 2;
            
                    
            mask = zeros(M,N);
            mask(x3:x4,y3:y4) = 1;
            
            fft_filter_holo = fft_holo .* mask;
            holo_filter_1 = fftshift(ifft2(fftshift(fft_filter_holo)));
        
            fft_holo_2 = fftshift(fft2(fftshift(holo(:,:,2)))); 
            fft_filter_holo_2 = fft_holo_2 .* mask;
            holo_filter_2 = fftshift(ifft2(fftshift(fft_filter_holo_2)));
        
        
            holo_filter(:,:,1) = holo_filter_1;
            holo_filter(:,:,2) = holo_filter_2;
        
            holo_FT(:,:,1) = fft_filter_holo;
            holo_FT(:,:,2) = fft_filter_holo_2;
        
            if visual == 'Yes'
                figure,imagesc(log(abs(fft_filter_holo).^2)),colormap(gray),title('FT Filter Hologram'),daspect([1 1 1]) 
                figure,imagesc((abs(holo_filter_1).^2)),colormap(gray),title('Filter Hologram'),daspect([1 1 1]) 
            end
        end


      
        function [cf] = cost_function(seed_maxPeak,lambda,dxy,M,N,m,n,k,fx_0,fy_0,holo_filtered)
            cf = 0;
            fx_max = seed_maxPeak(1,1);
            fy_max = seed_maxPeak(1,2);
            k = 2*pi/lambda;
            theta_x = asin((fx_0 - fx_max) * lambda / (M * dxy));
            theta_y = asin((fy_0 - fy_max) * lambda / (N * dxy));
            ref_wave = exp(1i * k * (sin(theta_x) * n * dxy + sin(theta_y) * m * dxy));
        
            holo_rec = holo_filtered .* ref_wave;
            phase = angle(holo_rec);
        
            phase_save = uint8(255 * mat2gray(phase));
            ib = imbinarize(phase_save, 0.1);
            cf = M*N - sum(ib(:));
        end


        function [cf] = costFunction_SIDHM(theta, FTHolo, fx_max,fy_max)
            cf = 0;
            [M,N] = size(FTHolo);
            
            [Dtemp] = functions_evaluation.demComp2SIDHM(theta,FTHolo);
            Dplus = Dtemp(:,:,1);
            Dplus = abs(Dplus).^0.1;

            maxL = abs(Dplus(fy_max(1),fx_max(1)));
            maxR = abs(Dplus(fy_max(2),fx_max(2)));
            
            %cf = 1/abs(maxR - maxL);
            cf =  maxL/abs(maxR + maxL);
        end

        function [cf] = costFunction_SIDHM_II(theta, FTHolo, fx_max,fy_max)
            cf = 0;
            [M,N] = size(FTHolo);
            
            [Dtemp] = functions_evaluation.demComp2SIDHM(theta,FTHolo);%demodulation in the Fourier domain
            Dplus = Dtemp(:,:,2);
            Dplus = abs(Dplus).^0.1;
            
            maxL = abs(Dplus(fy_max(1),fx_max(1)));
            maxR = abs(Dplus(fy_max(2),fx_max(2)));
            
            cf =  maxR/abs(maxR + maxL);
    
        end

        function [D] = demComp2SIDHM(theta, H)
            [X, Y, no] = size(H); 
            D = zeros(X,Y,no);
            M = 1/2*[exp(1i*theta(1)) exp(-1i*theta(1));exp(1i*theta(2)) exp(-1i*theta(2))];
            Minv = inv(M);
            
            D(:,:,1) = Minv(1,1).*H(:,:,1) + Minv(1,2).*H(:,:,2);
            D(:,:,2) = Minv(2,1).*H(:,:,1) + Minv(2,2).*H(:,:,2);
        end
    end
end


