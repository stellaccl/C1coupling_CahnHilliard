%phtSquarePlateFixed.m
%
% PHT Isogeometric analysis for Kirchoff plate problems.
%
% Rotation-free thin plates. Fully clamped rectangular plates.
%
%

addpath('./exampleData')
addpath('./utils')
close all
clear

tic;

numPatches = 1;
W = 100;
L = 100;
GIFTmesh = init2DGeometryGIFTMP('rectangle', L, W, numPatches);

dimBasis = zeros(1, numPatches);
p=3;
q=3;

PHUTelem = cell(numPatches, 1);
for i=1:numPatches
    [PHUTelem{i}, dimBasis(i)] = initPHTmesh(p,q);
    quadList{i} = 2:5;
end

% Material properties

E  = 1e7;
nu = 0.3;
t  = L/1000; % thickness

% Boundary condition

q0  = -1.;  % distributed force

clamped = 1;

% constitutive matrix

D  = E*t^3/(12*(1-nu^2));
C  = D*[1 nu 0;nu 1 0; 0 0 0.5*(1-nu)];

%elasticity matrix (plane stress)
Cmat=E/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];

%define the number of refinement steps
%numSteps =6;
numSteps=3;
for stepCount = 1:numSteps
    
    figure
    plotPHTMeshMP(PHUTelem, GIFTmesh)
    
    patchBoundaries = [];
    [PHUTelem, sizeBasis ] = zipConforming_IGAPack( PHUTelem, dimBasis, patchBoundaries, p, q); 
    [ stiff, rhs ] = assembleGalerkinPlate( PHUTelem, GIFTmesh, sizeBasis, p, q, C, q0);

    disp('Imposing boundary conditions...')
    %[ stiff, rhs, bcdof, bcval ] = imposeDirichletNeuPHMPPlate(stiff, rhs, PHUTelem, p, q);
    [ stiff, rhs, bcdof, bcval ] = imposeDirichletNeuPHMPPlate2(stiff, rhs, PHUTelem, p, q);
    % pause
    %toc
    disp('Solving the linear system...')
    condest(stiff)
    sol0 = stiff\rhs;

    minZ = PlotDispPlate( PHUTelem, GIFTmesh,  p, q, sol0);
    
    quadRef = ones(1,size(quadList{1},1));
    [quadList{1}, PHUTelem{1}, dimBasis] = refineMesh(quadRef, quadList{1}, PHUTelem{1}, p, q, dimBasis);
    
    wbar = minZ*D*1000/q0/W^4
    wext = 1.26532  % CCCC

end


