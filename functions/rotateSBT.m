function [rotatedSbtM] = rotateSBT(SbtM,es)
% ROTATESBT
% Last updated: 05/10/2020
%--------------------------------------------------------------------------
% Instructions for use
%--------------------------------------------------------------------------
%ROTATESBT for generalising self-interaction matrix
%   Turns matrix SbtM which was valid for body-frame coordinates "es"
%   into matrix rotatedSbtM which is valid in the lab frame.
%--------------------------------------------------------------------------

if abs(det(es)-1) > 1e-3
   warning('Body-frame matrix is no longer orthogonal.') 
end

N = round(size(SbtM,1)/3);
backrotate = kron(eye(N),es);
rotate = kron(eye(N),inv(es));

rotatedSbtM = backrotate*SbtM*rotate;

end

