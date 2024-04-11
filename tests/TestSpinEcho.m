% Test script for BlochSim_CK_3D
% Single channel spin echo simulation
% M Zhang

clear all;
dt = 1e-5;
nPoints = ceil(1e-3/dt);

RFexc = ones(nPoints, 1);
RF = 2 * ones(nPoints, 1);
g = 12* 1e-3 * ones(nPoints, 1);
 
x = 0;
y = 0;
z = 0.1;
b1 = 250 / (42.5e6);
M0 = [0; 0; 1];
T1 = 10^20;
T2 = 10^20;
df = 0;

%% Build sequence
nT = floor(0.007/1e-5);
totalRF = zeros(nT, 1);
totalg = zeros(nT, 3);

totalRF(1:nPoints, 1) = RFexc;
totalRF(3*nPoints+1:4*nPoints, 1)= RF;
totalg(nPoints+1:2*nPoints, 3) = g;
totalg(5*nPoints+1:6*nPoints, 3) = g;
% Exc-Crush-Refoc-Crush-

%% Simulate
nz = 21;
zRange = linspace(-0.001,0.001,nz);

[Mxy, Mz] = blochSim_CK_3D(totalRF, totalg, dt, 50*ones(1,1,nz), x,y,zRange, b1*ones(1,1,nz), ...
    saveall=true);

%% plot
figure
subplot(4,1,1)
plot(totalRF)
xticks([])
ylim([0 2.2])
title('a) RF amplitude / a.u.')
subplot(4,1,2)
plot(totalg)
xticks([])
ylim([0 0.013])
title('b) Gradient / T/m')
subplot(4,1,3)
plot(squeeze(angle(Mxy)))
title('c) Phase')
xticks([])
subplot(4,1,4)
plot(1e-3:1e-2:1e-2*nT,squeeze(abs(mean(Mxy,4))))
title('d) Mean M_{xy}')
xlabel('t / ms')
