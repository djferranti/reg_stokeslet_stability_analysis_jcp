function parElasticNonlinearAlgebraicAperiodic(dt, tFinal, epsilonFac, N, bending)

%remove the references to the mat/file/vidObj/ if you do not wish to save 
%the data/text/video outputs
if bending 
    matname = ['aperiodic' num2str(N) 'BendingAlgebraicEpsilon' num2str(epsilonFac) 'h.mat'];
    fileID = fopen(['aperiodicBendingAlgebraicL2_' num2str(N) '_Epsilon' num2str(epsilonFac) '.txt'],'w');
    vidObj = VideoWriter(['aperiodicBendingAlgebraic' num2str(N) '_Epsilon' num2str(epsilonFac) '.avi'], "Motion JPEG AVI");
    vidObj.Quality = 100;
    open(vidObj);
else 
    matname = ['aperiodic ' num2str(N) 'TensionOnlyAlgebraicEpsilon' num2str(epsilonFac) 'h.mat'];
    fileID = fopen(['aperiodicTensionAlgebraicL2_' num2str(N) '_Epsilon' num2str(epsilonFac) '.txt'],'w');
    vidObj = VideoWriter(['aperiodicTensionAlgebraic' num2str(N) '_Epsilon' num2str(epsilonFac) '.avi'], "Motion JPEG AVI");
    vidObj.Quality=100;
    open(vidObj);
end

%% vidparams
skip = 1;
%% geometry
%physical grid parameters
L = 1;  
h = L / N;

%% fluid/simulation parameters
epsilon = epsilonFac * h; %regularization
mu = 1; %viscosity

%% elastic membrane initialization
K = 75; %tensile rigidity
if bending
    Kb = 0.001 * K; %bending rigidity
else 
    Kb = 0;
end

A0 = 0.01;  %initial perturbed amplitude of membrane
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
% zPerturbation = A0 .* sin( 2 .* pi./ L  .* (Xq + Xs) );
Xs = Xs + sPerturbation;
Xz = Xz + zPerturbation;

l2UPrev = 10^3;

%% simulate the perturbed membrane
timeStepCounter = 0;
increaseCounter = 0;
for t = 0 : dt : tFinal %do the following at every time step

    %compute force/geometry data for membrane
    [Fq, Fs, Fz, dSmembrane] = ...
        computetensilebendingforcesaperiodic(K, Xq, ...
        Xs, Xz, M1, M2, dq, ds, Kb); 

    [UMembrane, VMembrane, WMembrane] = ...
        evaluatevelocityaperiodic(Fq, Fs, Fz, Xq, Xs, Xz, dSmembrane, ...
        epsilon, mu);

    Xq = Xq + UMembrane .* dt;
    Xs = Xs + VMembrane .* dt;
    Xz = Xz + WMembrane .* dt;

    l2U = sqrt ( sum ( sum ( (UMembrane.^2 + VMembrane .^2 + ...
        WMembrane .^2) ./ (M1 .* M2) ) ) );

    if (mod(timeStepCounter, skip) == 0)
        vidObj=plotmembrane(vidObj,Xq, Xs, Xz, L, l2U, t);
    end 

    if (l2U>l2UPrev)
        increaseCounter = increaseCounter + 1;
    else 
        increaseCounter = 0;
    end 

    if l2U > 100 
        break
    end

    l2UPrev = l2U; 

    fprintf(fileID, 't = %4.4f \n', t);
    fprintf(fileID, 'l2 U = %.7e \n', l2U);

    %iterate time step counter
    timeStepCounter = timeStepCounter + 1;
end 
save(matname)
fclose(fileID);
close(vidObj);
end
%% helper functions

function [UMembrane, VMembrane, WMembrane] = ...
    evaluatevelocityaperiodic(Fq, Fs, Fz, ...
    Xq, Xs, Xz,dSmembrane, epsilon, mu)

UMembrane = zeros(size(Xq)); VMembrane = zeros(size(Xs));
WMembrane = zeros(size(Xz));

numberForces = numel(Fq);
parfor k = 1 : numberForces
 %for k = 1 : numberForces 
    forceK = [Fq(k), Fs(k), Fz(k)]' .* dSmembrane(k);
    
    [Uk, Vk, Wk] = regularizedstokesletclassic(Xq, Xs, Xz, Xq(k), ...
    Xs(k), Xz(k), forceK, epsilon, mu);


    UMembrane = UMembrane + Uk;
    VMembrane = VMembrane + Vk;
    WMembrane = WMembrane + Wk;
end
end

function [Uij, Vij, Wij] = regularizedstokesletclassic(X, Y, Z, x0ij, ...
    y0ij, z0ij, forceK, reg, mu)

Rx = X - x0ij; Ry = Y - y0ij; Rz = Z - z0ij;
R = (Rx.^2 + Ry.^2 + Rz.^2 + reg.^2).^(1/2);
Rm1 = R.^(-1); Rm3 = R.^(-3);
fDotR = forceK(1) .* Rx + forceK(2) .* Ry + forceK(3) .* Rz;

H1 = Rm1 + reg.^2 .* Rm3; 
H2 = Rm3; 

Uij = H1 .* forceK(1) + H2 .* fDotR .* Rx; 
Vij = H1 .* forceK(2) + H2 .* fDotR .* Ry;
Wij = H1 .* forceK(3) + H2 .* fDotR .* Rz; 

Uij = Uij ./ ( 8 * pi * mu ); 
Vij = Vij ./ ( 8 * pi * mu ); 
Wij = Wij ./ ( 8 * pi * mu ); 

end 

function vidObj = plotmembrane(vidObj, Xq, Xs, Xz, L, l2U, t)
figure(1)
hold off

surf(Xq,Xs,Xz,'FaceColor','interp','LineStyle','none')
hold on

skp = 1;
plot3(Xq(1:skp:end),Xs(1:skp:end),Xz(1:skp:end), ...
    '.k','LineWidth',1,'MarkerSize', 5)
view(-28, 38)
xlabel('x','FontSize',12)
ylabel('y','FontSize',12)
zlabel('z','FontSize',12)
xlim([-0.1,L+0.1])
ylim([-0.1,L+0.1]);
zlim([-0.02, 0.02])
title([' t = ' num2str(t) ', $\ell_2(\mathbf{U})$ = ' num2str(l2U)], ...
'Interpreter','latex', 'FontSize',14) ;
currFrame=getframe(gcf);
writeVideo(vidObj,currFrame);
end