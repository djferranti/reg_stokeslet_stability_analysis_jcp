function parElasticNonlinearEwaldGauss(dt, tFinal, epsilonFac, N, bending)
%implements the doubly periodic elastic surface-fluid problem using either
%of the Gaussian blobs, phi^{G,1M} or phi^{G,3M}. 
%% geometry
NInterp = 2 * N;

if bending 
    matname = ['parEwald ' num2str(N) 'BendingGaussianEpsilon' num2str(epsilonFac) 'h.mat'];
    fileID = fopen(['bendingGaussianL2_' num2str(N) '_Epsilon' num2str(round(epsilonFac)) '.txt'],'w');
    vidObj = VideoWriter(['bendingGaussian' num2str(N) '_Epsilon' num2str(epsilonFac) '.avi'], "Motion JPEG AVI");
    vidObj.Quality=100;
    open(vidObj);
else 
    matname = ['parEwald ' num2str(N) 'TensionOnlyGaussianEpsilon' num2str(epsilonFac) 'h.mat'];
    fileID = fopen(['tensionOnlyGaussianL2_' num2str(N) '_Epsilon' num2str(round(epsilonFac)) '.txt'],'w');
    vidObj = VideoWriter(['tensionGaussian' num2str(N) '_Epsilon' num2str(epsilonFac) '.avi'], "Motion JPEG AVI");
    vidObj.Quality=100;
    open(vidObj);
end

%% geometry 
skip =1; %for video plotting 
%physical grid parameters
L = 1;  %length of periodic box
hCoarse = L / N;
h = L / (NInterp); %fine grid size
xj = 0 : h : L-h; yj = xj; %grid points

%fourier grid parameters
kx = [ 0 : 1 : N / 2 - 1, N / 2, -N / 2 + 1 : 1 : -1]; %1d wavenumbers for FFT
ky = kx;
[Kx, Ky] = meshgrid(kx, ky); %2d wavenumbers in meshgrid style
Kx = Kx .* 2 * pi / L; Ky = Ky .* 2 * pi / L; %scaled fourier wave numbers for cont FT
C = sqrt(Kx.^2 + Ky.^2); %used in Fourier Green's functions

%% fluid/simulation parameters
epsilon = epsilonFac * hCoarse; %regularization
xi = 6 * hCoarse; %ewald splitting parameter
rCutoff = 1.75*L; %evaluate local velocity only within rCutoff

mu = 1; %viscosity

%% elastic surface initialization
K = 75; %tensile rigidity
Kb = 0.01 * K; %bending rigidity 
A0 = 0.01;  %initial perturbed amplitude of surface
dq = L / N;
ds = L / N;
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
%physical grid in meshgrid format
[X,Y,Z] = meshgrid(xj, yj, zLevels); %grid points in 2d array meshgrid style

isperiodic = 1;
%% simulate the perturbed surface
timeStepCounter = 0;
increaseCounter = 0;
for t = 0 : dt : tFinal %do the following at every time step

    %compute force/geometry data for surface
    if bending
        [FqEwald, FsEwald, FzEwald, dSsurfaceEwald] = ...
            computetensilebendingforces(K, Xq, ...
            Xs, Xz, N, N, dq, ds, L, bending, Kb);
    else
        [FqEwald, FsEwald, FzEwald, dSsurfaceEwald] = ...
            computetensilebendingforces(K, Xq, ...
            Xs, Xz, N, N, dq, ds, L);
    end

    %uses the gauss blob with one moment condition, phi^{G,1M}
    [U, V, W, UsurfaceLocal, VsurfaceLocal, WsurfaceLocal] = ...
        evaluatevelocityewaldgauss1m(FqEwald, FsEwald, FzEwald, Xq, Xs, Xz, ...
        X, Y, Z, dSsurfaceEwald,L, Kx, Ky, C, epsilon, xi, rCutoff, mu, NInterp);

    %uses the gauss blob with three moment conditions, phi^{G,3M}
    % [U, V, W, UsurfaceLocal, VsurfaceLocal, WsurfaceLocal] = ...
    %     evaluatevelocityewaldgauss3m(FqEwald, FsEwald, FzEwald, Xq, Xs, Xz, ...
    %     X, Y, Z, dSsurfaceEwald,L, Kx, Ky, C, epsilon, xi, rCutoff, mu, NInterp);

    %for debugging interp
    [UsurfaceLong, VsurfaceLong, WsurfaceLong] = ...
        interpolatetosurface(xj, yj, zLevels, L, U, V, ...
        W, Xq, Xs, Xz, isperiodic);

    Usurface = UsurfaceLong + UsurfaceLocal;
    Vsurface = VsurfaceLong + VsurfaceLocal;
    Wsurface = WsurfaceLong + WsurfaceLocal;

    Xq = Xq + Usurface .* dt;
    Xs = Xs + Vsurface .* dt;
    Xz = Xz + Wsurface .* dt;

    l2UEwald = sqrt ( sum ( sum ( (Usurface.^2 + Vsurface .^2 + ...
        Wsurface .^2) ./ (N .* N) ) ) );

    fprintf(fileID, 't = %4.4f \n', t);
    fprintf(fileID, 'l2 U = %.7e \n', l2UEwald);
    
    if (l2UEwald>l2UEwaldPrev)
        increaseCounter = increaseCounter + 1;
    else 
        increaseCounter = 0;
    end

    if (mod(timeStepCounter, skip) == 0)
        vidObj=plotsurface(vidObj,Xq, Xs, Xz, L, l2UEwald, t);
    end 

    %if you want to stop once L2 U increased in 10 sequential steps,
    %uncomment this

    % if increaseCounter > 9
    %     disp('L2 U surface increased in 10 sequential steps')
    %     disp(' ')
    %     isstable = 0; 
    %     break
    % end 

    if l2UEwald > 20
        break
    end

    l2UEwaldPrev = l2UEwald;


    %make a new z-levels vector based on the current surface config
    zLevels = initializezlevels(Xz, numberZLevels, isplanar);
    [X,Y,Z] = meshgrid(xj, yj, zLevels); %grid points in 2d array meshgrid style
    % [XLocal,YLocal,ZLocal] = meshgrid(xjlocal, yjlocal, zLevels); %grid points in 2d array meshgrid style

    %iterate time step counter
    timeStepCounter = timeStepCounter + 1;
end 
save(matname)
fclose(fileID);
close(vidObj);

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
    zLevels = [zMin, zMin + minSpacing];
end
end

function [Usurface, Vsurface, Wsurface] = interpolatetosurface(xj, ...
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
    Usurface = interp2( Xj, Yj, U, Xq, Xs );
    Vsurface = interp2( Xj, Yj, V, Xq, Xs );
    Wsurface = interp2( Xj, Yj, W, Xq, Xs );
else
    Usurface = interp3( Xj, Yj, Zj, U, Xq, Xs, Xz, 'spline' );
    Vsurface = interp3( Xj, Yj, Zj, V, Xq, Xs, Xz, 'spline' );
    Wsurface = interp3( Xj, Yj, Zj, W, Xq, Xs, Xz, 'spline' );
end

end

function [U2dp, V2dp, W2dp, UsurfaceLocal, VsurfaceLocal, WsurfaceLocal] = ...
    evaluatevelocityewaldgauss1m(Fq, Fs, Fz, ...
    Xq, Xs, Xz, X, Y, Z, dSsurface, L, Kx, Ky, C, ...
    epsilon, xi, rCutoff, mu, NInterp)

U2dp = zeros(size(X)); V2dp = zeros(size(Y)); W2dp = zeros(size(Z));
UsurfaceLocal = zeros(size(Xq)); VsurfaceLocal = zeros(size(Xs));
WsurfaceLocal = zeros(size(Xz));

numberForces = numel(Fq);
%use for if you don't want parallelization
parfor k = 1 : numberForces
    forceK = [Fq(k), Fs(k), Fz(k)]' .* dSsurface(k);
    xyzK = [Xq(k), Xs(k), Xz(k)]';
    [U2dpK, V2dpK, W2dpK, UkLocal, VkLocal, WkLocal] = ...
        dpregularizedstokesletewaldgauss1m(C, Kx, Ky, ...
        X, Y, Z, Xq, Xs, Xz, L, xyzK, forceK, epsilon, xi, rCutoff, ...
        mu, NInterp);

    U2dp = U2dp + U2dpK;
    V2dp = V2dp + V2dpK;
    W2dp = W2dp + W2dpK;

    UsurfaceLocal = UsurfaceLocal + UkLocal;
    VsurfaceLocal = VsurfaceLocal + VkLocal;
    WsurfaceLocal = WsurfaceLocal + WkLocal;
end
end

function [U2dp, V2dp, W2dp, UsurfaceLocal, VsurfaceLocal, WsurfaceLocal] = ...
    evaluatevelocityewaldgauss3m(Fq, Fs, Fz, ...
    Xq, Xs, Xz, X, Y, Z, dSsurface, L, Kx, Ky, C, ...
    epsilon, xi, rCutoff, mu, NInterp)

U2dp = zeros(size(X)); V2dp = zeros(size(Y)); W2dp = zeros(size(Z));
UsurfaceLocal = zeros(size(Xq)); VsurfaceLocal = zeros(size(Xs));
WsurfaceLocal = zeros(size(Xz));

numberForces = numel(Fq);
%use for if you don't want parallelization
parfor k = 1 : numberForces
    forceK = [Fq(k), Fs(k), Fz(k)]' .* dSsurface(k);
    xyzK = [Xq(k), Xs(k), Xz(k)]';
    [U2dpK, V2dpK, W2dpK, UkLocal, VkLocal, WkLocal] = ...
        dpregularizedstokesletewaldgauss3m(C, Kx, Ky, ...
        X, Y, Z, Xq, Xs, Xz, L, xyzK, forceK, epsilon, xi, rCutoff, ...
        mu, NInterp);

    U2dp = U2dp + U2dpK;
    V2dp = V2dp + V2dpK;
    W2dp = W2dp + W2dpK;

    UsurfaceLocal = UsurfaceLocal + UkLocal;
    VsurfaceLocal = VsurfaceLocal + VkLocal;
    WsurfaceLocal = WsurfaceLocal + WkLocal;
end
end

function vidObj = plotsurface(vidObj, Xq, Xs, Xz, L, l2U, t)

figure(1)
hold off
Xq = [Xq, Xq + L]; Xq = repmat(Xq, 2, 1);
Xs = [Xs; Xs + L]; Xs = repmat(Xs, 1, 2); 
Xz = repmat(Xz,2,2);
surf(Xq,Xs,Xz,'FaceColor','interp','LineStyle','none')
hold on

skp = 1;
plot3(Xq(1:skp:end),Xs(1:skp:end),Xz(1:skp:end), ...
    '.k','LineWidth',1,'MarkerSize', 5)
view(-28, 38)
xlabel('x','FontSize',12)
ylabel('y','FontSize',12)
zlabel('z','FontSize',12)
xlim([-0.1,2*L+0.1])
ylim([-0.1,2*L+0.1]);
zlim([-0.02, 0.02])
title([' t = ' num2str(t) ', $\ell_2(\mathbf{U})$ = ' num2str(l2U)], ...
'Interpreter','latex', 'FontSize',14) ;
currFrame=getframe(gcf);
writeVideo(vidObj,currFrame);
end
end