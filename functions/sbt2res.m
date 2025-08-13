function [FullRes] = sbt2res(FullSBT,es,fs,epsil,psi,Nturns,c,NLegendre)
% SBT2RES  
% Last updated: 02/05/2020
%--------------------------------------------------------------------------
% Instructions for use
%--------------------------------------------------------------------------
% SBT2RES Go from SBT matrix to resistance matrix.
% INPUTS
%   FUllSBT     full SBT matrix in terms of Legendre polynomial modes
%   es          orientation of first filament
%   fs          orientation of second filament
%   epsil       slenderness parameter
%   psi         helical angle (set to 0 if straight)
%   Nturns      number of turns (set to 1 if straight)
%   c           chirality of the helix (defaults to +1, right-handed)
%   NLegendre   number of Legendre modes (defaults to 10)
%
% WORKING VARIABLES
%   rjk     function, centreline position of j'th filament, k'th Cartesian component
%   tjk     function, centreline tangent of j'th filament, k'th Cartesian component
%   u1      3-by-1 vector function, velocity of points along first filament
%   u2      3-by-1 vector function, velocity of points along second filament
%   f       3-by-N matrix, force distribution along filament;
%           first dim is for cartesian components,
%           second dim is for expanding in terms of Legendre polynomials
%   F       3-by-1 row vector, total force
%   T       3-by-1 row vector, total torque
%
%   Pn        N-by-1 cell array of Legendre polynomials
%   Pnvector  N-by-1 vector function of Legendre polynomials
%
% OUTPUTS
%   FullRes   6*NFilaments-by-3*NFilaments extended resistance matrix
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Legendre polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MLegendre = zeros(NLegendre);
MLegendre(1,1) = 1;
MLegendre(2,2) = 1;
for p=3:NLegendre
    MLegendre(p,:) = (2*p-3)/(p-1)*[0 MLegendre(p-1,1:NLegendre-1)] - (p-2)/(p-1)*MLegendre(p-2,:);
end
Pnvector = @(s) MLegendre*s.^transpose(0:NLegendre-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Commonly used quantities
%%%%%%%%%%%%%%%%%%%%%%%%%%%
piN  = pi*Nturns;
SIN  = sin(psi);
COS  = cos(psi);
R  = SIN/piN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define centreline, tangent, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First filament (relative to centre of roation)
r11 = @(s) R*cos(piN*s)*es(1,1)+c*R*sin(piN*s)*es(1,2)+COS*s*es(1,3);
r12 = @(s) R*cos(piN*s)*es(2,1)+c*R*sin(piN*s)*es(2,2)+COS*s*es(2,3);
r13 = @(s) R*cos(piN*s)*es(3,1)+c*R*sin(piN*s)*es(3,2)+COS*s*es(3,3);

t11 = @(s) -SIN*sin(piN*s)*es(1,1)+c*SIN*cos(piN*s)*es(1,2)+COS*s.^0*es(1,3);
t12 = @(s) -SIN*sin(piN*s)*es(2,1)+c*SIN*cos(piN*s)*es(2,2)+COS*s.^0*es(2,3);
t13 = @(s) -SIN*sin(piN*s)*es(3,1)+c*SIN*cos(piN*s)*es(3,2)+COS*s.^0*es(3,3);

% Second filament (relative to centre of roation)
r21 = @(s) R*cos(piN*s)*fs(1,1)+c*R*sin(piN*s)*fs(1,2)+COS*s*fs(1,3);
r22 = @(s) R*cos(piN*s)*fs(2,1)+c*R*sin(piN*s)*fs(2,2)+COS*s*fs(2,3);
r23 = @(s) R*cos(piN*s)*fs(3,1)+c*R*sin(piN*s)*fs(3,2)+COS*s*fs(3,3);

t21 = @(s) -SIN*sin(piN*s)*fs(1,1)+c*SIN*cos(piN*s)*fs(1,2)+COS*s.^0*fs(1,3);
t22 = @(s) -SIN*sin(piN*s)*fs(2,1)+c*SIN*cos(piN*s)*fs(2,2)+COS*s.^0*fs(2,3);
t23 = @(s) -SIN*sin(piN*s)*fs(3,1)+c*SIN*cos(piN*s)*fs(3,2)+COS*s.^0*fs(3,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute resistance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise output
FullRes = zeros(12);

% tic
for k=1:6
    % First let helix 1 move at this velocity
    
    % Loop over linear velocities and angular velocities in principal
    % directions x,y,z
    U = double(k==[1; 2; 3]);
    Omega = double(k==[4; 5; 6]);
    u1 = @(s) U + cross(Omega,[r11(s); r12(s); r13(s)],1);
    
    % Calculate velocity components
    % First helix
    fun = @(s) 8*pi*kron(Pnvector(s),u1(s));
    a1 = integral(fun,-1,1,'ArrayValued',true);
    a1 = reshape(a1,3*NLegendre,1);
    % Second helix
    a2 = zeros(3*NLegendre,1);
    % Stick together
    a = [a1; a2];
    
    % Calculate force components
    f = FullSBT\a;
    f1 = f(1:3*NLegendre);
    f2 = f((3*NLegendre+1):end);
    f1 = reshape(f1,3,NLegendre);
    f2 = reshape(f2,3,NLegendre);
    
    % Compute total force
    F1 = 2*f1(:,1);
    F2 = 2*f2(:,1);
    
    % Compute total torque
    T1 = zeros(3,1);
    T2 = zeros(3,1);
    for n = 1:NLegendre
        fun = @(s) [r11(s); r12(s); r13(s)] .* (MLegendre(n,:)*s.^transpose(0:NLegendre-1));
        T1 = T1 + cross(integral(fun,-1,1,'ArrayValued',true),f1(:,n),1);
        
        fun = @(s) [r21(s); r22(s); r23(s)] .* (MLegendre(n,:)*s.^transpose(0:NLegendre-1));
        T2 = T2 + cross(integral(fun,-1,1,'ArrayValued',true),f2(:,n),1);
    end
    
    if k>=4
        % Add contribution from rotlet singularities
        rot = @(s) 4*pi*epsil^2*(1-s^2)*...
            dot(Omega,[t11(s); t12(s); t13(s)],1)*...
            [t11(s); t12(s); t13(s)];
        T1 = T1 + integral(rot,-1,1,'ArrayValued',true);
    end
    
    % Enter force and torque into resistance matrix
    FullRes( 1:3, k) = F1;
    FullRes( 4:6, k) = T1;
    FullRes( 7:9, k) = F2;
    FullRes(10:12,k) = T2;
      
    % Now let helix 2 move at this velocity
    
    % Loop over linear velocities and angular velocities in principal
    % directions x,y,z
    U = double(k==[1; 2; 3]);
    Omega = double(k==[4; 5; 6]);
    u2 = @(s) U + cross(Omega,[r21(s); r22(s); r23(s)],1);
    % N.B. Careful to rotate second helix about its own axis
    
    % Calculate velocity components
    % First helix
    a1 = zeros(3*NLegendre,1);
    % Second helix
    fun = @(s) 8*pi*kron(Pnvector(s),u2(s));
    a2 = integral(fun,-1,1,'ArrayValued',true);
    a2 = reshape(a2,3*NLegendre,1);
    % Stick together
    a = [a1; a2];
    
    % Calculate force components
    f = FullSBT\a;
    f1 = f(1:3*NLegendre);
    f2 = f((3*NLegendre+1):end);
    f1 = reshape(f1,3,NLegendre);
    f2 = reshape(f2,3,NLegendre);
    
    % Compute total force
    F1 = 2*f1(:,1);
    F2 = 2*f2(:,1);
    
    % Compute total torque
    T1 = zeros(3,1);
    T2 = zeros(3,1);
    for n = 1:NLegendre
        fun = @(s) [r11(s); r12(s); r13(s)] .* (MLegendre(n,:)*s.^transpose(0:NLegendre-1));
        T1 = T1 + cross(integral(fun,-1,1,'ArrayValued',true),f1(:,n),1);
        
        fun = @(s) [r21(s); r22(s); r23(s)] .* (MLegendre(n,:)*s.^transpose(0:NLegendre-1));
        T2 = T2 + cross(integral(fun,-1,1,'ArrayValued',true),f2(:,n),1);
    end
    
    if k>=4
        % Add contribution from rotlet singularities
        rot = @(s) 4*pi*epsil^2*(1-s^2)*...
            dot(Omega,[t21(s); t22(s); t23(s)],1)*...
            [t21(s); t22(s); t23(s)];
        T2 = T2 + integral(rot,-1,1,'ArrayValued',true);
    end
    
    % Enter force and torque into resistance matrix
    FullRes( 1:3, k+6) = F1;
    FullRes( 4:6, k+6) = T1;
    FullRes( 7:9, k+6) = F2;
    FullRes(10:12,k+6) = T2;
end

end

