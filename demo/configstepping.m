% CONFIGSTEPPING v1pt0                                        13-AUG-2025 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GNU General Public License v3.0
% Author: Maria Tătulea-Codrean
% Reference: M. Tătulea-Codrean & E. Lauga. (2021) 
% Asymptotic theory of hydrodynamic interactions between slender filaments. 
% Physical Review Fluids, 6: 074103.
%--------------------------------------------------------------------------
%
% This routine computes the hydrodynamic forces and torques on two helical 
% filaments rotating in parallel in a viscous fluid, with a prescribed 
% angular velocity (configstepping = stepping through a prescribed sequence 
% of configurations at a constant rate). The configuration sequence
% can be adapted using the variables for position (x)
% and orientation (es) of the two filaments.
%
% Code structure:
% 1) Lines 34-59: Define physical parameters and sequence of configurations.
%
% 2) Lines 61-72: If needed, modify numerical parameters. 
%
% 3) Lines 74-172: Hydrodynamic computations based on Slender Body Theory. 
% This is the most computationally intensive step. It can take a few minutes 
% to solve the hydrodynamics problem for each configuration (repat for the
% number of configurations in the requested sequence). The code automatically 
% saves progress along the way. You can stop and restart the simulation 
% at any time, with minimal losses.
%
% 4) Lines 174-212: Interpret simulation results and generate plots.
%--------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-defined parameters [SI units]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'configstepping_sampledata.mat'; % change for each new simulation

% Filament geometry
r = 10*1e-3;    % helix radius [m]
p = 40*1e-3;    % helix pitch [m]
re = 2.5*1e-3;  % cross-section radius [m]
h = 100*1e-3;   % axial length [m]
c = -1;         % chirality (-1 for left-handed, +1 for right-handed)

% Physical properties
mu   = 1;       % dynamic viscosity [Pa s]
freq = 0.1;     % frequency [Hz]

% Configuration
d    = 30*1e-3;     % inter-axial distance [m]
dphi = 0;           % phase difference between filaments

% Decode inputs
[l, Omega, Fscale, Tscale] = generateDimensionalParameters(h,p,r,freq,mu);
[Nturns,psi,epsil,R,L,D] = generateDimensionlessParameters(h,p,r,re,l,d);

% Define configuration sequence
[x,es,~,~,phi1] = defineConfigurationSequence('dphi', dphi, 'distance', D, 'phi1', linspace(0,2*pi,11));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical accuracy parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N.B. increasing the accuracy (by decreasing the tolerance values or
% increasing the number of Legendre polynomials) will increase computation time

% Accuracy for computing double integrals
rtolr = 1e-5; % MATLAB default is 1e-6
atolr = 1e-6; % MATLAB default is 1e-10

% Legendre poly truncation level
NLegendre = 15; % recommended value: NLegendre > 3*Nturns

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfile(filename)
    % Save workspace data (inc. status)
    status = 'notstarted';
    save(filename)
    disp(['saving ' filename])
else
    % Read status
    load(filename,'status')
end

switch true
    case strcmp(status,'notstarted')
        % Simulation is about to start
        disp('initializing variables')

        % Initialise variables
        Nfilaments = 2; % number of filaments
        FullRes = zeros(Nfilaments*6,Nfilaments*6,size(x,4));
        start_kk = 1;
        save(filename,'Nfilaments','FullRes','-append')

        % Ready for hydrodynamic computations
        disp('starting simulation')

    case strcmp(status,'inprogress')
        % Simulation is in progress
        disp(['loading ' filename])

        % Load counter (kk) and output variable (FullRes)
        load(filename,'kk','FullRes')
        start_kk = kk;

        % Ready for hydrodynamic computations
        disp(['continuing simulation from counter ' num2str(start_kk) '/' num2str(size(x,4))])

    case strcmp(status,'done')
        % Simulation is finished
        disp('simulation completed')

        % Ready for analysing results
        disp('starting data analysis')
end

if ~strcmp(status,'done')
    % Look for previously computed ResM and SbtM
    if strcmp(status,'inprogress')
        load(filename,'ResM','SbtM')
    end

    % If necessary, compute ResM and SbtM for the first time
    if (~exist('ResM','var') || ~exist('SbtM','var'))
        % Intrinsic resistance/SBT matrix for isolated helix
        [ResM, SbtM] = sbtself(epsil,psi,Nturns,c,NLegendre,rtolr,atolr);
        save(filename,'ResM','SbtM','-append')
    end

    starttime = cputime;
    for kk = start_kk:size(x,4)
        % Simulation has started
        status = 'inprogress';
        save(filename,'status','kk','-append')

        % Read current configuration
        x0 = x(:,:,:,kk);
        es0 = es(:,:,:,kk);

        % Rotate SBT matrix into configuration of each helix
        SBT11 = rotateSBT(SbtM,es0(:,:,1));
        SBT22 = rotateSBT(SbtM,es0(:,:,2));

        % Compute cross-interaction matrix
        SBT12 = sbtcross(x0(:,:,1),x0(:,:,2),es0(:,:,1),es0(:,:,2),...
            epsil,psi,Nturns,c,NLegendre,rtolr,atolr);

        % Put SBT matrix together
        FullSBT = [SBT11, SBT12; SBT12', SBT22];

        % Integrate SBT matrix to get Res(istance) matrix
        FullRes(:,:,kk) = sbt2res(FullSBT,es0(:,:,1),es0(:,:,2),epsil,psi,Nturns,c,NLegendre);

        % Logging time
        disp([filename ' progress:' num2str(kk) '/' num2str(size(x,4)) ...
            ', running time: ' num2str((cputime-starttime)/60) ' min'])

        % Saving progress
        save(filename,'FullRes','-append')
    end

    % Simulation completed
    status = 'done';
    save(filename,'status','FullRes','-append')
    disp(['saving ' filename])

    disp('simulation completed')
    disp('starting data analysis')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: In this demo, the control variable is phi1, while d and dphi are 
% input parameters, but the configuration space could be defined by any
% combination of vectors of these three variables. Adjust plot as needed.
load(filename,'FullRes','phi1','d','dphi')

% Normal force & torque
Fz = Fscale*(FullRes(3,6,:)+FullRes(3,12,:)); % normal force on first filament, due to unit angular velocity rotation of *both* filaments
Tz = Tscale*(FullRes(6,6,:)+FullRes(6,12,:)); % axial torque on first filament, due to unit angular velocity rotation of *both* filaments

% Plot
mycols = lines(2);

figure('Position',[499,388,756,258])
subplot(1,2,1)
convf = 1e3; % unit converxion from N to mN
plot(phi1(:),convf*Fz(:),'--o','Color',mycols(1,:),'MarkerSize',8,'MarkerFaceColor',mycols(1,:))
xlabel('\phi_1')
ylabel('normal force, F_z [mN]')
box on
grid on
xlim([0 2*pi])
xticks([0 1 2]*pi)
xticklabels({'0','\pi','2\pi'})

subplot(1,2,2)
convf = 1e5; % unit converxion from N*m to mN*cm
plot(phi1(:),convf*Tz(:),'--o','Color',mycols(2,:),'MarkerSize',8,'MarkerFaceColor',mycols(2,:))
xlabel('\phi_1')
ylabel('normal force, T_z [mN cm]')
box on
grid on
xlim([0 2*pi])
xticks([0 1 2]*pi)
xticklabels({'0','\pi','2\pi'})

sgtitle(['\Delta\phi = ' num2str(dphi) ', d = ' num2str(d*1e3) ' mm'])