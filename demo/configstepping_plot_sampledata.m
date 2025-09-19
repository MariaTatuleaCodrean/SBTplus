% Note: In the sample data, the control variable is phi1, while d and dphi are 
% input parameters, but the configuration space could be defined by any
% combination of vectors of these three variables. Adjust plot as needed.
load('configstepping_sampledata.mat','FullRes','phi1','d','dphi','Fscale','Tscale')

% Normal force & torque
Fz = Fscale*FullRes(3,6,:);
Tz = Tscale*FullRes(6,6,:);

% Plot
mycols = lines(2);

figure('Position',[499,388,756,258])
subplot(1,2,1)
convf = 1e3;  % conversation from N to mN
plot(phi1(:),convf*Fz(:),'--o','Color',mycols(1,:),'MarkerSize',8,'MarkerFaceColor',mycols(1,:))
xlabel('\phi_1')
ylabel('normal force, F_z [mN]')
box on
grid on
xlim([0 2*pi])
xticks([0 1 2]*pi)
xticklabels({'0','\pi','2\pi'})

subplot(1,2,2)
convf = 1e5;  % conversation from N*m to mN*cm
plot(phi1(:),convf*Tz(:),'--o','Color',mycols(2,:),'MarkerSize',8,'MarkerFaceColor',mycols(2,:))
xlabel('\phi_1')
ylabel('normal force, T_z [mN cm]')
box on
grid on
xlim([0 2*pi])
xticks([0 1 2]*pi)
xticklabels({'0','\pi','2\pi'})

sgtitle(['\Delta\phi = ' num2str(dphi) ', d = ' num2str(d*1e3) ' mm'])
