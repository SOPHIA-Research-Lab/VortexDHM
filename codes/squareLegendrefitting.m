%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:-->        Square Legendre polynomials                        
%                                                                                                               
%                                                                                  
% Description:-->  This algorithm allows to calculate the first 15th square 
%                  Legendre polynomials.              
%                                                                                  
% Authors:-->      Rene Restrepo(2), Karina Ortega-Sanchez(1,2), 
%                  Carlos Trujillo(2), Ana Doblas(3) and Alfonso Padilla (1)
%                  (1) Polytechnic University of Tulancingo
%                  (2) EAFIT Univeristy SOPHIA Research Group (Applied Optics) 
%                  (3) Department of Electrical and Computer Engineering, University 
%                      of Massachusetts Dartmouth
%
% Inputs:-->       (1) Polynomial order
%                  (2) x and y coordinates
% Outputs:-->      (1) Square Legendre polynomials
%
% Email:-->        rrestre6@eafit.edu.co                                         
% Date:-->         12/17/2024                                                      
% version 1.0 (2024)                                                               
% Notes-->       
function polynomials = squareLegendrefitting(order, x, y)
all_polynomials = {1;
    x;
    y;
    ((3.*(x.^2))-1)/2;
    x.*y;
    ((3.*(y.^2))-1)/2;
    (x.*((5.*(x.^2))-3))/2;
    (y.*((3.*(x.^2))-1))/2;
    (x.*((3.*(y.^2))-1))/2;
    (y.*((5.*(y.^2))-3))/2;
    ((35.*(x.^4))-(30.*(x.^2))+3)/8;
    (x.*y.*((5.*(x.^2))-3))/2;
    (((3.*(y.^2))-1).*((3.*(x.^2))-1))/4;
    (x.*y.*((5.*(y.^2))-3))/2;
    ((35.*(y.^4))-(30.*(y.^2))+3)/8;};

polynomials = zeros(size(x));
for i = 1:length(order)
    polynomials(:,:,i) = all_polynomials{order(i)};
end
end