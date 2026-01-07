%% doubly periodic examples
%example with Gaussian type blob
dt = 7.2e-3; %time step 
tFinal = 5;  %simulate up to tFinal 
epsilonFac = 4; %use epsilon = epsilonFac*h where h=1/N
N = 32;  %grid points per direction
bending = 1; %use 1 for bending, 0 for surface with no bending

parElasticNonlinearEwaldGauss(dt, tFinal, epsilonFac, N, bending);

%example with algebraic blob
% dt = 3.7e-3; %time step 
% tFinal = 5;  %simulate up to tFinal 
% epsilonFac = (3/4*sqrt(pi))^(1/3)*2; %use epsilon = epsilonFac*h where h=1/N
% N = 32;  %grid points per direction
% bending = 1; %use 1 for bending, 0 for surface with no bending
% 
% parElasticNonlinearAlgebraic(dt, tFinal, epsilonFac, N, bending);

%% aperiodic example
% t = 3e-3; %time step 
% tFinal = 5;  %simulate up to tFinal 
% epsilonFac = (3/4*sqrt(pi))^(1/3)*2; %use epsilon = epsilonFac*h where h=1/N
% N = 32;  %grid points per direction
% bending = 1;
% 
% parElasticNonlinearAlgebraicAperiodic(dt, tFinal, epsilonFac, N, bending)
