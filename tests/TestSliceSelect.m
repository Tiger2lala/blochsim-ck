% Test script for BlochSim_CK_3D
% Single channel 90 degree slice selection pulse
% M Zhang

clear all;
addpath("..")

%% Parameters
dt = 1e-5; % s
nPoints = ceil(1e-3/dt);
 
RF = 1 * ones(nPoints, 1); % V
g = 40* 1e-3 * ones(nPoints, 1); % T/m
 
x = 0;
y = 0;
z = 0.1; % m
b1 = 250 / (42.5e6); % T/V, 1ms for 90 deg
M0 = [0; 0; 1];
T1 = 10^20; % s
T2 = 10^20;
df = 0; % Hz

%% Build sequence
nT = floor(1e-3/1e-5);
totalRF = zeros(nT, 1);
totalg = zeros(nT, 3);

totalRF(1:nPoints, 1)= RF;
totalg(1:nPoints, 3) = g;

%% Sim
nz = 21;
zRange = linspace(-0.001,0.001,nz);

[Mxy, Mz] = blochSim_CK_3D(totalRF, totalg, dt, zeros(1,1,nz), x,y,zRange, b1*ones(1,1,nz), ...
    saveall=true);

%%
figure
subplot(3,1,1)
plot(totalRF); ylabel("RF")
subplot(3,1,2)
plot(totalg); ylabel("Gz")
subplot(3,1,3)
plot(zRange, squeeze(abs(Mz(end,1,1,:)))); ylabel("Mz")

figure
subplot(2, 1, 1); hold on;
plot(abs(squeeze(Mxy)))
title("Mxy evolution of all voxels")

subplot(2, 1, 2)
plot((squeeze(Mz))); title("Mz evolution of all voxels")