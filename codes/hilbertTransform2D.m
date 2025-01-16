function cuadrature = hilbertTransform2D(c,HilbertOrEnergyOperator)

% Hilbert Transform 2D usually namaed Spiral Phase Transform
% cuadrature = hilbertTransform2D(c,1) calculates the two dimentional Hilbert transform.
% Given a matrix c=b*cos(psi), cuadrature is equal to cuadrature = i*exp(i*beta)*sin(psi) 
% where beta is the direction of the fringe map in every position.  

% Energy Operator relates to the two dimentional Hilbert Tranform.
% cuadrature = hilbertTransform2D(c,0)
% Given a matrix c=b*cos(psi), the result would be            
% energyOperator=-b*exp(i*beta)*sin(psi), where beta es the direccion of the firnge
% map.


% Kieran G. Larkin, Donald J. Bone, and Michael A. Oldfield,
% "Natural demodulation of two-dimensional fringe patterns. I.
% General background of the spiral phase quadrature transform,"
% J. Opt. Soc. Am. A 18, 1862-1870 (2001)

% Larkin, Kieran G. Uniform estimation of orientation using
% local and nonlocal 2-D Energy operators. 2005. Optical Society of
% America.



[NR, NC]=size(c);
[u,v]=meshgrid(1:NC, 1:NR);
u0=floor(NC/2)+1;
v0=floor(NR/2)+1;

u=u-u0;
v=v-v0;

H=(u+1i*v)./abs(u+1i*v);
H(v0, u0)=0; 

C=fft2(c);

if HilbertOrEnergyOperator
    CH=C.*ifftshift(H);
else
    CH=C.*ifftshift(1i.*H);
end

cuadrature = conj(ifft2(CH));
end
