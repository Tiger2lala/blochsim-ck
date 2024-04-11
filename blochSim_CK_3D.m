
function [Mxy, Mz] = blochSim_CK_3D(RF, g, dt, df, x, y, z, b1, options)
% Native matlab version
% Based on Pauly (1991) DOI:10.1109/42.75611
% M Zhang 2022
%
% Relaxation NOT yet included.
%
% Input:
%     'RF'    [nT, nCha]  RF-trajectory [V]
%     'g'     [nT, 3]     Gradient trajectory [T/m]
%     'dt'    [1]         Time step [s]
%     'df'    [ny,nx,nz]  Off resonance frequency [Hz]
%     'x'     [nx]        values of x [m]
%     'y'     [ny]        values of y [m]
%     'z'     [nz]        values of z [m]
%     'b1'    [ny,nx,nz,nCha] [T/V]
%     'M0'    [3,]        Initial magnetization
%     'T1'    [1]         [s]
%     'T2'    [1]         [s]
%     'saveall' 0 or 1
% Output:
%     'Mxy'   [nT,ny,nx,nz] Complex Mxy
%     'Mz'    [nT,ny,nx,nz] Real Mz



arguments
    RF                              % (nT, nCha)  RF-trajectory [V]
    g                               % (nT, 3)     Gradient trajectory [T/m]
    dt      double                  % (1, 1)      Time step [s]
    df                              % (ny,nx,nz)  Off resonance frequency [Hz]
    x                               % (1, nx)     position values of x [m]
    y
    z
    b1                              % (ny,nx,nz,nCha) [T/V]
    options.M0  (3,1)   = [0;0;1]   % Initial magnetization
    options.T1  double = 100        % [s]
    options.T2  double = 100        % [s]
    options.saveall logical = false % Save all steps
end

gamma = 267 * 10^6; % 1H [rad s-1 T-1]
gammabar = gamma/(2*pi); % [Hz T-1]
R1 = 1/options.T1; R2 = 1/options.T2;

%% Size and meshgrid
[nT, nCha] = size(RF);
if ~(size(b1,4)==nCha)
    error("No. of Channels in RF and b1 must match")
end

if ~(all(size(b1,[2 1 3])==[numel(x) numel(y) numel(z)]))
    error("Spatial locations in b1 must match x y z")
end

if ~(all(size(df,[2 1 3])==[numel(x) numel(y) numel(z)]))
    error("Spatial locations in df must match x y z")
end

[X,Y,Z] = meshgrid(x, y, z); % becomes [ny, nx, nz]

%% simulation
if options.saveall==1
    saveLen = nT;
else
    saveLen = 1;
end

finalM = complex(zeros(3,numel(X)));
A = ones(1,1,numel(X));
B = zeros(1,1,numel(X));
 
M0 = options.M0;
M = [M0(1)+1j*M0(2); M0(1)-1j*M0(2); M0(3)];
Mxy = complex(zeros(saveLen, numel(X))); 
Mz = complex(zeros(saveLen, numel(X)));

for iDx = 1 : nT

    RF_T = reshape(reshape(b1, [], nCha)*RF(iDx,:).', numel(y), numel(x), numel(z));
    Bx = real(RF_T(:,:,:)); % [ny, nx, nz] in T
    By = imag(RF_T(:,:,:));
    Bz = df/gammabar + X * g(iDx,1) + Y* g(iDx,2) + Z * g(iDx, 3); % in T
 
    Btot = cat(4, Bx, By, Bz); % [ny, nx, nz, 3]
    phi = -gamma * dt .* (sqrt( vecnorm(Btot, 2, 4).^2 ) + 1e-50); % rad
    n = (gamma * dt ./ abs(phi)) .* Btot;
 
    alpha = reshape(cos(phi/2) - 1i* n(:,:,:,3) .* sin(phi/2),1,1,[]); % [1,1,nVox]
    beta = reshape(-1i * (n(:,:,:,1) + 1i * n(:,:,:,2)) .* sin(phi/2), 1,1 ,[]);
 
    Anew = alpha.*A-conj(beta).*B;
    Bnew = beta.*A+conj(alpha).*B;
%     Q_all(iDx, :, :, :) = [alpha, -conj(beta); beta, conj(alpha)];
%     Q_step = [alpha, -conj(beta); beta, conj(alpha)]; % [2,2,nVox]
%  
%     Q = pagemtimes(Q_step,Q); % [2,1,nVox]
%  
%     alpha = Q(1, 1, :); % [1,1,nVox]
%     beta = Q(2, 1, :);
    A = Anew;
    B = Bnew;

    if M(1)==0 % only need 2 elements!
        finalM(1,:) = squeeze(2*conj(A).*B); %[1, nVox]
        finalM(3,:) = squeeze(A.*conj(A)-B.*conj(B)); 
    else
        convert_mat = [conj(A) .^2, -B.^2, 2 * conj(A) .* B; ...
            -(conj(B).^2), A.^2, 2 * A .* conj(B); ...
            -conj(A).*conj(B), -A.*B, A.*conj(A)-B.*conj(B)]; % [3,3,nVox]
    
        finalM = squeeze(pagemtimes(convert_mat, M)); % [3,nVox]
    end


    if options.saveall==1
        Mxy(iDx,:) = finalM(1,:); % [nT, nVox]
        Mz(iDx,:) = finalM(3,:);
    elseif iDx==nT
        Mxy(1,:) = finalM(1,:); % [1, nVox]
        Mz(1,:) = finalM(3,:);
    end

end

Mxy = reshape(Mxy, [saveLen size(X)]); 
Mz = reshape(Mz, [saveLen size(X)]);

