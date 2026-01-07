function isstable = parElasticNonlinearNewAlgebraic(dt, tFinal, epsilonFac, N, bending)
isstable=1;
%% geometry
%N = 32; %grid size in physical/fourier space
NInterp = 2 * N;

if bending 
    matname = ['par ' num2str(N) 'BendingNewAlgebraicEpsilon' num2str(round(epsilonFac)) 'h.mat'];
    fileID = fopen(['bendingNewAlgL2_' num2str(N) '_Epsilon' num2str(round(epsilonFac)) '.txt'],'w');
else 
    matname = ['par ' num2str(N) 'TensionOnlyAlgebraicEpsilon' num2str(round(epsilonFac)) 'h.mat'];
    fileID = fopen(['tensionOnlyNewAlgL2_' num2str(N) '_Epsilon' num2str(round(epsilonFac)) '.txt'],'w');
end
% matname='parEwald32GaussianEpsilon2h.mat';
%physical grid parameters
L = 1;  %length of periodic box
hCoarse = L / N;
h = L / (NInterp); %fine grid size
% hFine = L / NInterp; %mesh size
xj = 0 : h : L-h; yj = xj; %grid points

%fourier grid parameters
kx = [ 0 : 1 : N / 2 - 1, N / 2, -N / 2 + 1 : 1 : -1]; %1d wavenumbers for FFT
ky = kx;
[Kx, Ky] = meshgrid(kx, ky); %2d wavenumbers in meshgrid style
Kx = Kx .* 2 * pi / L; Ky = Ky .* 2 * pi / L; %scaled fourier wave numbers for cont FT
C = sqrt(Kx.^2 + Ky.^2); %used in Fourier Green's functions

%% fluid/simulation parameters
epsilon = epsilonFac * hCoarse; %regularization
% xi = 8 * hCoarse; %ewald splitting parameter
rCutoff = L; %evaluate local velocity only within rCutoff

mu = 1; %viscosity

skp=1;
tPoints = 0 : dt : tFinal;

%% elastic membrane initialization
K = 75; %tensile rigidity
% bending = 1; %membrane has bending rigidity
Kb = 0.01 * K; %bending rigidity 
A0 = 0.01;  %initial perturbed amplitude of membrane
% A0 = 0.005;
M1 = N;
M2 = N;
dq = L / M1;
ds = L / M2;
qj = 0 : dq : L-dq; sj = qj; %grid points
[Xq, Xs] = meshgrid(qj, sj);
Xz = zeros ( size (Xq) );

sPerturbation = A0 .* sin( 4 .* pi ./ L .* Xq );
zPerturbation = A0 .* sin( 2 .* pi./ L  .* (Xq + Xs) ) + ...
    A0 .* cos (4 .* pi./L .* Xq);
Xs = Xs + sPerturbation;
Xz = Xz + zPerturbation;

l2UEwaldPrev = 10^3;

numberZLevels = 21;
isplanar = 0;
zLevels = initializezlevels(Xz, numberZLevels, isplanar);
zMinInitial = min(Xz, [], 'all');
zMaxInitial = max(Xz, [], 'all');
%physical grid in meshgrid format
[X,Y,Z] = meshgrid(xj, yj, zLevels); %grid points in 2d array meshgrid style

isperiodic = 1;
%% simulate the perturbed membrane
timeStepCounter = 1;
increaseCounter = 0;
for t = 0 : dt : tFinal %do the following at every time step

    %compute force/geometry data for membrane
    [FqEwald, FsEwald, FzEwald, dSmembraneEwald] = ...
        computetensilebendingforces(K, Xq, ...
        Xs, Xz, M1, M2, dq, ds, L, bending, Kb);

    [U, V, W] = ...
        evaluatevelocity(FqEwald, FsEwald, FzEwald, Xq, Xs, Xz, ...
        X, Y, Z, dSmembraneEwald,L, Kx, Ky, C, epsilon, mu, NInterp);
    % evaluatevelocity(Fq, Fs, Fz, ...
    % Xq, Xs, Xz, X, Y, Z, dSmembrane, L, Kx, Ky, C, ...
    % epsilon, mu, NInterp)

    [UMembrane, VMembrane, WMembrane] = ...
        interpolatetomembrane(xj, yj, zLevels, L, U, V, ...
        W, Xq, Xs, Xz, isperiodic);

    Xq = Xq + UMembrane .* dt;
    Xs = Xs + VMembrane .* dt;
    Xz = Xz + WMembrane .* dt;

    l2UEwald = sqrt ( sum ( sum ( (UMembrane.^2 + VMembrane .^2 + ...
        WMembrane .^2) ./ (M1 .* M2) ) ) );
    %disp(['l2 = ' num2str(l2UEwald)]);
    if (l2UEwald>l2UEwaldPrev)
        %disp('l2 U increased')
        increaseCounter = increaseCounter + 1;
    else 
        increaseCounter = 0;
    end

    if increaseCounter > 9
        disp('L2 U membrane increased in 10 sequential steps')
        isstable = 0; 
        break
    end
    l2UEwaldPrev = l2UEwald;

    fprintf(fileID, 't = %4.4f \n', t);
    fprintf(fileID, 'l2 U = %.7e \n', l2UEwald);


    %make a new z-levels vector based on the current membrane config
    zLevels = initializezlevels(Xz, numberZLevels, isplanar);
    [X,Y,Z] = meshgrid(xj, yj, zLevels); %grid points in 2d array meshgrid style

    %iterate time step counter
    timeStepCounter = timeStepCounter + 1;
end 
save(matname)
fclose(fileID);
end
%% helper functions

function zLevels = initializezlevels(Xz, numberZLevels, isplanar)
zMin = min(Xz, [], 'all');
zMax = max(Xz, [], 'all');
minSpacing = 1e-4;

if ~isplanar && (zMax-zMin)>minSpacing
    minZLevels = max( ceil( (zMax - zMin) /  minSpacing ), 5);
    numberZLevels = min(numberZLevels, minZLevels);
    zLevels = [linspace(zMin, zMax, numberZLevels), zMax + minSpacing];
else
    zLevels = [zMin, zMax];
end
end

function [Umembrane, Vmembrane, Wmembrane] = interpolatetomembrane(xj, ...
    yj, zLevels, L, U, V, W, Xq, Xs, Xz, isperiodic)

%the following takes care of the periodicity in the interpolation by adding
%L in the x,y directions.
%otherwise, a point may be closer to the end of the periodic boundary than
%than the endpoint of our grid and will not be interpolated.

if isperiodic
    [Xj, Yj, Zj] = meshgrid( [xj, L], [yj, L], zLevels);
    U = cat( 2, U, U(:, 1, :) );
    U = cat( 1, U, U(1, :, :) );
    V = cat( 2, V, V(:, 1, :) );
    V = cat( 1, V, V(1, :, :) );
    W = cat( 2, W, W(:, 1, :) );
    W = cat( 1, W, W(1, :, :) );
else
    [Xj, Yj, Zj] = meshgrid( xj, yj, zLevels );
end

Xq = mod( Xq, L);
Xs = mod( Xs, L);

if zLevels == 0
    Umembrane = interp2( Xj, Yj, U, Xq, Xs );
    Vmembrane = interp2( Xj, Yj, V, Xq, Xs );
    Wmembrane = interp2( Xj, Yj, W, Xq, Xs );
else
    Umembrane = interp3( Xj, Yj, Zj, U, Xq, Xs, Xz, 'spline' );
    Vmembrane = interp3( Xj, Yj, Zj, V, Xq, Xs, Xz, 'spline' );
    Wmembrane = interp3( Xj, Yj, Zj, W, Xq, Xs, Xz, 'spline' );
end

end

function [U2dp, V2dp, W2dp] = ...
    evaluatevelocity(Fq, Fs, Fz, ...
    Xq, Xs, Xz, X, Y, Z, dSmembrane, L, Kx, Ky, C, ...
    epsilon, mu, NInterp)

U2dp = zeros(size(X)); V2dp = zeros(size(Y)); W2dp = zeros(size(Z));

numberForces = numel(Fq);
% tic
for k = 1 : numberForces
    forceK = [Fq(k), Fs(k), Fz(k)]' .* dSmembrane(k);
    xyzK = [Xq(k), Xs(k), Xz(k)]';
    [U2dpK, V2dpK, W2dpK] = ...
        dpregularizedstokesletnewalgebraic(C, Kx, Ky, ...
        X, Y, Z, L, xyzK, forceK, epsilon, ...
        mu, NInterp);

    U2dp = U2dp + U2dpK;
    V2dp = V2dp + V2dpK;
    W2dp = W2dp + W2dpK;

end
% toc
end


