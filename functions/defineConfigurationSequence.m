function [x,es,DVEC,DPHI,PHI1,PHI2] = defineConfigurationSequence(varargin)
%DEFINECONFIGURATIONSEQUENCE Define position and orientation arrays
% Last updated: 13/08/2025
% Usage:
%   [x,es] = defineConfigurationSequence('dphi', pi/2, 'distance', 2,'phi1', [0 pi/2 pi 3*pi/2 2*pi]);
%   [x,es] = defineConfigurationSequence('distance', 2, 'phi1', 0, 'dphi', linspace(0,pi,6))
%
% Inputs 
%   dphi        phase difference between filaments [rad]
%   distance    dimensionless inter-axial distance
%   phi1        rotation angle of first filament about z-axis
%
% Outputs
%   x       position of filaments (standard form, see below)
%   es      orientation of filaments (standard form, see below)
%
% *for each configuration in the sequence*
%   DVEC    vector of dimensionless inter-axial distance 
%   DPHI    vector of phase difference 
%   PHI1    vector of first filament phase angle
%   PHI2    vector of second filament phase angle
%
%--------------------------------------------------------------------------
% Standard form (for x, es)
% dim = 1   Cartesian coordinates
% dim = 2   multiple columns for matrices
% dim = 3   filament index
% dim = 4   time/configuration index
%--------------------------------------------------------------------------

% Create input parser
p = inputParser;

% Define default values
defaultDistance = 2;
defaultDphi = 0;
defaultPhi1 = linspace(0,2*pi,13);

% Add parameters with validation
addParameter(p, 'distance', defaultDistance, @(x) isa(x,'double'));
addParameter(p, 'dphi', defaultDphi, @(x) isa(x,'double'));
addParameter(p, 'phi1', defaultPhi1, @(x) isa(x,'double'));

% Parse the inputs
parse(p, varargin{:});

% Extract parsed results
d = p.Results.distance;
dphi = p.Results.dphi;
phi1 = p.Results.phi1;

% Configuration grid
[D,DPHI,PHI1] = ndgrid(d,dphi,phi1);

% Phase angle of second filament
PHI2 = PHI1 + DPHI;

% Vectorize
DVEC = D(:);
DPHI = DPHI(:);
PHI1 = PHI1(:);
PHI2 = PHI2(:);
nconfig = length(D);

% Define output variables
x  = nan(3,1,2,nconfig);
es = nan(3,3,2,nconfig);

% Define position
for kk=1:nconfig
    % First filament (at origin)
    x(:,:,1,:)  = 0;
    
    % Second filament
    x(1,:,2,kk) = D(kk);
    x(2,:,2,kk) = 0;
    x(3,:,2,kk) = 0;
end

% Define orientation
for kk=1:nconfig
    % Orientation matrix
    % - each column of es corresponds to unit basis vector (e1, e2, e3) in
    % a body-fixed frame of reference

    % First filament
    es(:,:,1,kk) = [cos(PHI1(kk)), -cos(PHI1(kk)-pi/2), 0; cos(PHI1(kk)-pi/2), cos(PHI1(kk)), 0; 0, 0, 1];

    % Second filament
    es(:,:,2,kk) = [cos(PHI2(kk)), -cos(PHI2(kk)-pi/2), 0; cos(PHI2(kk)-pi/2), cos(PHI2(kk)), 0; 0, 0, 1];
end

end