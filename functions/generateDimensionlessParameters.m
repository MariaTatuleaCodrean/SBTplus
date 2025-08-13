function [Nturns,psi,epsil,R,L,D] = generateDimensionlessParameters(h,p,r,re,l,d)
% GENERATEDIMENSIONLESSPARAMETERS     
% Last updated: 13/08/2025
%--------------------------------------------------------------------------
% Instructions for use
%--------------------------------------------------------------------------
% Inputs
%   h       axial length [m]
%   p       helix pitch [m]
%   r       helix radius [m]
%   re      cross-sectional radius [m]
%   l       contour length of filament [m]
%   d       inter-axial distance [m]
%
% Outputs
%   Nturns  number of helical turns
%   psi     pitch angle of helix (relative to axis)
%   epsil   filament aspect ratio
%   R       dimensionless helix amplitude
%   L       dimensionless contour length
%   D       dimensionless inter-axial distance
%--------------------------------------------------------------------------

Nturns  = h/p;
psi     = atan(2*pi*r/p);
epsil   = 2*re/l;
R       = 2*r/l; % equivalent to sin(psi)/(pi*Nturns)
L       = 2*l/l;
D       = 2*d/l;

end