function [ResM,SbtM,fmom] = sbtself(epsil,psi,Nturns,c,NLegendre,rtolr,atolr)
% SBTSELF
% Last updated: 20/10/2020
%--------------------------------------------------------------------------
% History
%--------------------------------------------------------------------------
%   SBT v3.0 updated 20/10/2020 (calculate force moment F1k as well)
%   SBT v2.0 finalised 15/07/2019 (note: no actual change took place on 29
% June 2020, I just commented and then uncommented some lines)
%   SBT v1.0 created 08/07-12/07/2019 (discrete filaments, not working)
%   SBT v0.3 updated 03/07/2019
%   SBT v0.2 updated 26/03/2019
%   SBT v0.1 created 25/03/2019
%
%--------------------------------------------------------------------------
% Instructions for use
%--------------------------------------------------------------------------
% INPUTS
%   epsil      slenderness parameter
%   psi        helical angle (set to 0 if straight)
%   Nturns     number of turns (set to 1 if straight)
%
% OPTIONAL
%   c          chirality of the helix (defaults to +1, right-handed)
%   NLegendre  number of Legendre modes (defaults to 5)
%   rtolr
%   atolr
%
% WORKING VARIABLES
%   u       3-by-1 vector function, velocity of points along filament
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
%   ResM    6-by-6 resistance matrix
%   SbtM    3*NLegendre-by-3*NLegendre SBT matrix
%   fmom    6-by-1 force moment, useful for my asymptotic theory for hrydro
%   interactions between slender filaments
%
% Restrictions
%   1) This code only works with a straight or helical filament
%--------------------------------------------------------------------------

%%%%%%%%%%%%
% Preamble
%%%%%%%%%%%%
% Check inputs
if nargin<3
    error('Give at least curve and slenderness parameter.');
end
if nargin<4 || isempty(c)
    c = -1;
end
if nargin<5 || isempty(NLegendre)
    NLegendre = 5;
end
if nargin<6 || isempty(rtolr)
    rtolr = 1e-3;
end
if nargin<7 || isempty(atolr)
    atolr = 1e-4;
end

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Legendre polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MLegendre = zeros(NLegendre);
MLegendre(1,1) = 1;
MLegendre(2,2) = 1;

for k=3:NLegendre
    MLegendre(k,:) = (2*k-3)/(k-1)*[0 MLegendre(k-1,1:NLegendre-1)] - (k-2)/(k-1)*MLegendre(k-2,:);
end

Pnvector = @(s) MLegendre*s.^transpose(0:NLegendre-1);

%%%%%%%%%%%%%%%%
% Define shape
%%%%%%%%%%%%%%%%
piN  = pi*Nturns;
piN2 = piN/2;
SIN  = sin(psi);
SIN2 = SIN^2;
COS  = cos(psi);
COS2 = COS^2;
SINCOS = SIN*COS;
R  = SIN/piN;
R2 = R^2;
RR = 2*R;
RR2 = RR^2;
r1 = @(s) R*cos(piN*s);
r2 = @(s) c*R*sin(piN*s);
r3 = @(s) COS*s;
t1 = @(s) -SIN*sin(piN*s);
t2= @(s) c*SIN*cos(piN*s);
t3 = @(s) COS*s.^0;
titj = @(s) [SIN^2*sin(piN*s)^2,...
    -SIN*sin(piN*s)*c*SIN*cos(piN*s),...
    -SIN*sin(piN*s)*COS;...
    -SIN*sin(piN*s)*c*SIN*cos(piN*s),...
    SIN^2*cos(piN*s)^2,...
    c*SIN*cos(piN*s)*COS;...
    -SIN*sin(piN*s)*COS,...
    c*SIN*cos(piN*s)*COS,...
    COS^2*s^0];

%%%%%%%%%%%%%%%%%%%%%
% Compute SBT matrix
%%%%%%%%%%%%%%%%%%%%%
% Initialise output
SbtM = zeros(3*NLegendre);

for n=0:NLegendre-1
    % Determine eigenvalue for n
    En = 2*sum(1./(1:n));
    
    % Determine M^3_m,n factor (since it only depends on n)
    M3 = -En+log(4/epsil^2)-3;
    
    % Determine M^2_m,n (since it only depends on n)
    M2 = 2*(M3+4)/(2*n+1);
    
    % Add local contribution
    loc = 3*n+(1:3);
    SbtM(loc,loc) = SbtM(loc,loc) + M2*eye(3);
    
    % Calculate contribution from M^3_m,n,i,j
    fun = @(s) kron(Pnvector(s),titj(s))*(MLegendre(n+1,:)*s.^transpose(0:NLegendre-1));
    SbtM(:,loc) = SbtM(:,loc) + M3*integral(fun,-1,1,'ArrayValued',true);
    
    for m=0:NLegendre-1
        % Keep singularities on the integration boundaries
        sdlim = @(sd) sd;
        
        % 11
        ii = 1; jj = 1;
        
        % Define the integrand
        fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:NLegendre-1),size(s)).*...
            reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:NLegendre-1),size(sd)).*...
            (1./sqrt(RR2*sin(piN2*(s-sd)).^2+COS2*(s-sd).^2)+...
            (R*cos(piN*s)-R*cos(piN*sd)).^2./...
            (RR2*sin(piN2*(s-sd)).^2+COS2*(s-sd).^2).^(3/2)-...
            (1+SIN2*sin(piN*s).^2)./abs(s-sd));
        
        % Calculate contribution from M^1_m,n,i,j
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + integral2(fun,-1,1,-1,sdlim,'RelTol',rtolr,'AbsTol',atolr);
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + integral2(fun,-1,1,sdlim,+1,'RelTol',rtolr,'AbsTol',atolr);
        
        % 22
        ii = 2; jj = 2;
        
        % Define the integrand
        fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:NLegendre-1),size(s)).*...
            reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:NLegendre-1),size(sd)).*...
            (1./sqrt(RR2*sin(piN2*(s-sd)).^2+COS2*(s-sd).^2)+...
            (R*sin(piN*s)-R*sin(piN*sd)).^2./...
            (RR2*sin(piN2*(s-sd)).^2+COS2*(s-sd).^2).^(3/2)-...
            (1+SIN2*cos(piN*s).^2)./abs(s-sd));
        
        % Calculate contribution from M^1_m,n,i,j
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + integral2(fun,-1,1,-1,sdlim,'RelTol',rtolr,'AbsTol',atolr);
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + integral2(fun,-1,1,sdlim,+1,'RelTol',rtolr,'AbsTol',atolr);
        
        % 33
        ii = 3; jj = 3;
        
        % Define the integrand
        fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:NLegendre-1),size(s)).*...
            reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:NLegendre-1),size(sd)).*...
            (1./sqrt(RR2*sin(piN2*(s-sd)).^2+COS2*(s-sd).^2)+...
            (COS*s-COS*sd).^2./...
            (RR2*sin(piN2*(s-sd)).^2+COS2*(s-sd).^2).^(3/2)-...
            (1+COS2)./abs(s-sd));
        
        % Calculate contribution from M^1_m,n,i,j
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + integral2(fun,-1,1,-1,sdlim,'RelTol',rtolr,'AbsTol',atolr);
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + integral2(fun,-1,1,sdlim,+1,'RelTol',rtolr,'AbsTol',atolr);
        
        % 12 and 21
        ii = 1; jj = 2;
        
        % Define the integrand
        fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:NLegendre-1),size(s)).*...
            reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:NLegendre-1),size(sd)).*...
            (c*R2*(cos(piN*s)-cos(piN*sd)).*(sin(piN*s)-sin(piN*sd))./...
            (RR2*sin(piN2*(s-sd)).^2+COS2*(s-sd).^2).^(3/2)-...
            (0-c*SIN2*sin(piN*s).*cos(piN*s))./abs(s-sd));
        
        % Calculate contribution from M^1_m,n,i,j
        I = integral2(fun,-1,1,-1,sdlim,'RelTol',rtolr,'AbsTol',atolr);
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + I;
        SbtM(3*m+jj,3*n+ii) = SbtM(3*m+jj,3*n+ii) + I;
        
        I = integral2(fun,-1,1,sdlim,+1,'RelTol',rtolr,'AbsTol',atolr);
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + I;
        SbtM(3*m+jj,3*n+ii) = SbtM(3*m+jj,3*n+ii) + I;
        
        % 13 and 31
        ii = 1; jj = 3;
        
        % Define the integrand
        fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:NLegendre-1),size(s)).*...
            reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:NLegendre-1),size(sd)).*...
            (R*COS*(cos(piN*s)-cos(piN*sd)).*(s-sd)./...
            (RR2*sin(piN2*(s-sd)).^2+COS2*(s-sd).^2).^(3/2)-...
            (0-SINCOS*sin(piN*s))./abs(s-sd));
        
        % Calculate contribution from M^1_m,n,i,j
        I = integral2(fun,-1,1,-1,sdlim,'RelTol',rtolr,'AbsTol',atolr);
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + I;
        SbtM(3*m+jj,3*n+ii) = SbtM(3*m+jj,3*n+ii) + I;
        
        I = integral2(fun,-1,1,sdlim,+1,'RelTol',rtolr,'AbsTol',atolr);
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + I;
        SbtM(3*m+jj,3*n+ii) = SbtM(3*m+jj,3*n+ii) + I;
        
        % 23 and 32
        ii = 2; jj = 3;
        
        % Define the integrand
        fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:NLegendre-1),size(s)).*...
            reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:NLegendre-1),size(sd)).*...
            (c*R*COS*(sin(piN*s)-sin(piN*sd)).*(s-sd)./...
            (RR2*sin(piN2*(s-sd)).^2+COS2*(s-sd).^2).^(3/2)-...
            (0+c*SINCOS*cos(piN*s))./abs(s-sd));
        
        % Calculate contribution from M^1_m,n,i,j
        I = integral2(fun,-1,1,-1,sdlim,'RelTol',rtolr,'AbsTol',atolr);
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + I;
        SbtM(3*m+jj,3*n+ii) = SbtM(3*m+jj,3*n+ii) + I;
        
        I = integral2(fun,-1,1,sdlim,+1,'RelTol',rtolr,'AbsTol',atolr);
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + I;
        SbtM(3*m+jj,3*n+ii) = SbtM(3*m+jj,3*n+ii) + I;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute resistance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise output
ResM = zeros(6);
fmom = zeros(1,6);

for k=1:6
    % Loop over linear velocities and angular velocities in principal
    % directions x,y,z
    U = double(k==[1; 2; 3]);
    Omega = double(k==[4; 5; 6]);
    u = @(s) U + cross(Omega,[r1(s); r2(s); r3(s)],1);
    
    % Calculate velocity components
    fun = @(s) 8*pi*kron(Pnvector(s),u(s));
    a = integral(fun,-1,1,'ArrayValued',true);
    a = reshape(a,3*NLegendre,1);
    
    % Calculate force components
    f = SbtM\a;
    f = reshape(f,3,NLegendre);
    
    % Compute total force
    F = 2*f(:,1);
    
    % Compute total torque
    T = zeros(3,1);
    for n = 1:NLegendre
        fun = @(s) [r1(s); r2(s); r3(s)] .* (MLegendre(n,:)*s.^transpose(0:NLegendre-1));
        T = T + cross(integral(fun,-1,1,'ArrayValued',true),f(:,n),1);
    end
    
    if k>=4
        % Add contribution from rotlet singularities
        rot = @(s) 4*pi*epsil^2*(1-s^2)*...
            dot(Omega,[t1(s); t2(s); t3(s)],1)*...
            [t1(s); t2(s); t3(s)];
        T = T + integral(rot,-1,1,'ArrayValued',true);
    end
    
    ResM(1:3,k) = F;
    ResM(4:6,k) = T;
    
    % Compute force moment F1k
    for n = 1:NLegendre
        fun = @(s) [r1(s); r2(s); r3(s)] .* (MLegendre(n,:)*s.^transpose(0:NLegendre-1));
        fmom(k) = fmom(k) + dot(integral(fun,-1,1,'ArrayValued',true),f(:,n),1);
        
        fun = @(s) r1(s) * (MLegendre(n,:)*s.^transpose(0:NLegendre-1));
        fmom(k) = fmom(k) - 3*(integral(fun,-1,1,'ArrayValued',true))*f(1,n);
    end
    
end
now = toc;
disp(['Compute SBT self-resistance matrix - finished:' num2str(now) 's'])
end

