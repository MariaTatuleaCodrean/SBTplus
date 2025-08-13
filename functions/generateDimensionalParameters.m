function [l, Omega, Fscale, Tscale] = generateDimensionalParameters(h,p,r,freq,mu)
% GENERATEDIMENSIONALPARAMETERS     
% Last updated: 13/08/2025
%--------------------------------------------------------------------------
% Instructions for use
%--------------------------------------------------------------------------
% Inputs
%   h       axial length [m]
%   p       helix pitch [m]
%   r       helix radius [m]
%   freq    frequency [Hz]
%   mu      dynamic viscosity [Pa s]
%
% Outputs
%   l       contour length of the filament [m]
%   Omega   angular velocity [rad /s]
%   Fscale  force scale (multiply dimensionless SBTplus outputs by this)
%   Tscale  torque scale (multiply dimensionless SBTplus outputs by this)
%--------------------------------------------------------------------------

l = h/p*sqrt((2*pi*r)^2+p^2);   
Omega = 2*pi*freq; 
Fscale = Omega*mu*(l/2)^2;
Tscale = Omega*mu*(l/2)^3;

end