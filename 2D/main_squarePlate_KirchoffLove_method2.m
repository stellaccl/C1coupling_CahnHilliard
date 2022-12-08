restoredefaultpath
addpaths

close all
clear all

tic;

numPatches = 2;
W = 100;
L = 100;

GIFTmesh = init2DGeometryGIFTMP('rectangle', L, W, numPatches);

p=2;
q=2;

numElemU=2;
numElemV=2;

PHUTelem = cell(numPatches, 1);
dimBasis = zeros(1, numPatches);
quadList = cell(numPatches, 1);
for indexPatch=1:numPatches
    [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmesh_quadratic( p,q,numElemU,numElemV);
    quadList{indexPatch} =1:4;
end

patchBoundaries = {1,2,2,4};

%define the number of refinement steps
numSteps =5;
for stepCount = 1:numSteps
    
%     figure
%     plotPHTMeshMP(PHUTelem, GIFTmesh)
%     pause
    
    [PHUTelem,  sizeBasis ] = zipConformingNew( PHUTelem, dimBasis, patchBoundaries, p, q);
    
    disp('assign sol Index')
    [PHUTelem,solIndexCount] = assignSolIndex2D_quadratic(PHUTelem,patchBoundaries,sizeBasis,p,q);
    
%     figure
%     plotPHTMesh_solIndex( PHUTelem,GIFTmesh,p)
%     title('sol Index')
%     pause
%     
%     figure
%     plotPHTMesh_nodesGlobal( PHUTelem,GIFTmesh,p)
%     title('nodes Global')
%     pause
    
    disp('compute matrix m')
    [PHUTelem,m] =zipConformingC1_quadratic(PHUTelem,GIFTmesh,patchBoundaries,solIndexCount,p,q);
    
    disp('apply boundary condition')
    [coefSol] = zipConformingC1_BC_squarePlateKL(PHUTelem,m,p,q);
    
    disp('assign new nodes global and modify c')
    [PHUTelem,sizeBasis,type2Basis] = zipConformingC1_quaddratic_modifyC(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,coefSol,solIndexCount,p,q);
    
    disp('localising type 2 basis function')
    [PHUTelem] = localizeType2Basis2(PHUTelem,type2Basis);
    
    %     for i=1:sizeBasis
    %
    %         testPlotBasisPhys_GIFTmesh3(PHUTelem,GIFTmesh,p,q,i)
    %
    %         pause
    %         close all
    %     end
    
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
    
    disp('assembleGalerkin')
    [ stiff, rhs ] = assembleGalerkinPlate_c1( PHUTelem, GIFTmesh, sizeBasis, p, q, C, q0);
    
    disp('Imposing boundary conditions...')
    [ stiff, rhs, bcdof, bcval ] = imposeDirichletNeuPHMPPlate2Patch_c1(stiff, rhs, PHUTelem, p, q, numPatches,type2Basis);
    
    %toc
    disp('Solving the linear system...')
%     condest(stiff)
    sol0 = stiff\rhs;
    minZ =PlotDispSquarePlate_c1( PHUTelem, GIFTmesh,  p, q, sol0);
%     minZ = PlotDispPlate_c1( PHUTelem, GIFTmesh,  p, q, sol0)
    wbar = minZ*D*1000/q0/W^4
    wext = 1.26532
    pause
    %%%%%%%%%%%%%%%%%%%% refinement %%%%%%%%%%%%%%%%%%%%
    if stepCount<numSteps
        numElemU=numElemU*2;
        numElemV=numElemV*2;
        for indexPatch=1:numPatches
            [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmesh_quadratic( p,q,numElemU,numElemV);
            numElem=length(PHUTelem{1});
            quadList{indexPatch} =1:numElem;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end