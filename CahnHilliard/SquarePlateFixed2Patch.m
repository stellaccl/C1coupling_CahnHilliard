%phtSquarePlateFixed2Patch.m
%
% PHT Isogeometric analysis for Kirchoff plate problems.
%
% Rotation-free thin plates. Fully clamped rectangular plates.
% 2 patches
%

addpaths
addpath ./exampleData

close all
clear all

tic;

numPatches = 2;
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

%initialize the degree elevated mesh
% PHUTelemDE = cell(numPatches, 1);
% for i=1:numPatches
%     [PHUTelemDE{i}, dimBasisDE(i)] = initPHTmesh(p+1,q+1);
%     quadList{i} = 2:5;
% end


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
numSteps=1;

%do 1 uniform refinement
for indexPatch = 1:numPatches
    quadRef{indexPatch} = ones(1,size(quadList{indexPatch},1));
    [quadList{indexPatch}, PHUTelem{indexPatch}, dimBasis(indexPatch)] = refineMesh(quadRef{indexPatch}, quadList{indexPatch}, PHUTelem{indexPatch}, p, q, dimBasis(indexPatch));
end



for stepCount = 1:numSteps
    
    patchBoundaries = {1,2,2,4};
    %    pause
   
    [PHUTelem, sizeBasis ] = zipConforming_IGAPack( PHUTelem, dimBasis, patchBoundaries, p, q);
    
    figure
    plotPHTMeshMP(PHUTelem, GIFTmesh)

    %     globalType2basisNodes=1:sizeBasis;
%      testPlotBasisPhys_GIFTmesh3_degElev( PHUTelem,GIFTmesh,p,q,sizeBasis,globalType2basisNodes)
    %     pause
    %
    % %patchBoundaries = {1,2,2,4};
     patchBoundaries=[1,2,2,4];
    
    %rotation
    centerUpIndexPatch1=4;
    centerDownIndexPatch1=2;
    centerUpIndexPatch2=3;
    centerDownIndexPatch2=1;
    
    [ GIFTmesh ] = rotation( GIFTmesh,centerUpIndexPatch1,centerDownIndexPatch1,centerUpIndexPatch2,centerDownIndexPatch2);
    disp('zip conforming c1')
    [PHUTelem,sizeBasis,numType2Basis,type2Basis] =zipConforming_c1_arbitraryDegree_clamped(PHUTelem,GIFTmesh,sizeBasis, patchBoundaries, p, q );
    
    disp('Assembling LHS and RHS') 
    [ stiff, rhs ] = assembleGalerkinPlate_c1( PHUTelem, GIFTmesh, sizeBasis, p, q, C, q0);
    
    disp('Imposing boundary conditions...')
    [ stiff, rhs, bcdof, bcval ] = imposeDirichletNeuPHMPPlate2Patch_c1(stiff, rhs, PHUTelem, p, q, numPatches,type2Basis);
    % pause
    %toc
    disp('Solving the linear system...')
    condest(stiff)
    sol0 = stiff\rhs;
    
    minZ = PlotDispPlate_c1( PHUTelem, GIFTmesh,  p, q, sol0);
    disp('before refinement')
    pause
    %uniform refinement
    for indexPatch = 1:numPatches
        quadRef{indexPatch} = ones(1,size(quadList{indexPatch},1));
        [quadList{indexPatch}, PHUTelem{indexPatch}, dimBasis(indexPatch)] = refineMesh(quadRef{indexPatch}, quadList{indexPatch}, PHUTelem{indexPatch}, p, q, dimBasis(indexPatch));
    end
    wbar = minZ*D*1000/q0/W^4
    wext = 1.26532  % CCCC
end


