%phtSquarePlateFixed2Patch.m
% PHT Isogeometric analysis for Kirchoff plate problems.
% Rotation-free thin plates. Fully clamped rectangular plates.
% 2 patches
restoredefaultpath
addpaths

close all
clear all

tic;

numPatches = 2;
W = 100;
L = 100;
GIFTmesh = init2DGeometryGIFTMP('rectangle', L, W, numPatches);

p=3;
q=3;
numElemU=1;
numElemV=1;
PHUTelem = cell(numPatches, 1);
dimBasis = zeros(1, numPatches);
for indexPatch=1:numPatches
    [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmesh(p,q,numElemU,numElemV);
    quadList{indexPatch} = 2:5;
end

patchBoundaries = {1,2,2,4};

%define the number of refinement steps
numSteps =3;
for stepCount = 1:numSteps
    
    [PHUTelem,  sizeBasis ] = zipConformingNew( PHUTelem, dimBasis, patchBoundaries, p, q);
    
    [PHUTelem ,solIndexCount,patchInfo] = assignSolIndex2D( PHUTelem,patchBoundaries,p,q);
    
    figure
    plotPHTMeshMP(PHUTelem, GIFTmesh)
    %
    %     figure
    %     plotPHTMesh_nodesGlobal(PHUTelem, GIFTmesh,p)
    %     title('nodes Global')
    %
    %     figure
    %     plotPHTMesh_solIndex( PHUTelem,GIFTmesh,p)
    %     title('sol Index')
    %     pause
    
    disp('computing matrix m')
    [m] = zipConformingC1_computeM(PHUTelem,GIFTmesh,patchBoundaries,solIndexCount,p,q);
    
    disp('solving for coef sol')
    [coefSol] = zipConformingC1_boundaryCondition_SquarePlateClamped(PHUTelem,m,p,q );
    
    disp('modifying c')
    [ PHUTelem,sizeBasis,type2Basis] = zipConformingC1_modifyC(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,coefSol,p,q);
    
    disp('localising type 2 basis function')
    [PHUTelem] = localizeType2Basis2(PHUTelem,type2Basis);
    
    %%%%%% analysis part ===========================================
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
    %Cmat=E/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];

    disp('assembleGalerkin')
    [ stiff, rhs ] = assembleGalerkinPlate_c1( PHUTelem, GIFTmesh, sizeBasis, p, q, C, q0);
    
    disp('Imposing boundary conditions...')
    [ stiff, rhs, bcdof, bcval ] = imposeDirichletNeuPHMPPlate2Patch_c1(stiff, rhs, PHUTelem, p, q, numPatches,type2Basis);
    
    %toc
    disp('Solving the linear system...')
    %condest(stiff)
    sol0 = stiff\rhs;
    minZ =PlotDispSquarePlate_c1( PHUTelem, GIFTmesh,  p, q, sol0);
    %minZ = PlotDispPlate_c1( PHUTelem, GIFTmesh,  p, q, sol0)
    
    %uniform refinement
    for patchIndex=1:numPatches
        quadRef{patchIndex} = 1:size(quadList{patchIndex},1);
        [quadList{patchIndex}, PHUTelem{patchIndex}, dimBasis(patchIndex)] = refineMesh(quadRef{patchIndex}, quadList{patchIndex}, PHUTelem{patchIndex}, p, q, dimBasis(patchIndex));
    end
    
    
    format  long
    stepCount
    size_sol=size(sol0)
    wbar = minZ*D*1000/q0/W^4
    wext = 1.26532  % CCCC
end


