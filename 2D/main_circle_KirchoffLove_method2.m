close all
clear all
restoredefaultpath
addpaths
tic
numPatches = 5;
L=1;
W=1;

GIFTmesh = init2DGeometryGIFTMP('5PatchCircle',W,L,numPatches);

p=4;
q=4;

numElemU=4;
numElemV=4;
PHUTelem = cell(numPatches, 1);
dimBasis = zeros(1, numPatches);
quadList = cell(numPatches, 1);
if p==2
    for indexPatch=1:numPatches
        [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmesh_quadratic( p,q,numElemU,numElemV);
        %quadList{indexPatch} =1:4;
    end
else
    for i=1:numPatches
        [PHUTelem{i}, dimBasis(i)] = initPHTmeshGen( p,q,numElemU,numElemV);
        %quadList{i} = 1:4;
    end
end


patchBoundaries = {1,2,2,2;...
    2,3,4,2;...
    [1,3],4,[4,4],[2,4];...
    [1,2,3,4],5,[3,3,3,3],[2,3,4,1]};

numSteps =2;
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
%     %
%     figure
%     plotPHTMesh_nodesGlobal( PHUTelem,GIFTmesh,p)
%     title('nodes Global')
%     pause
    
    disp('compute matrix m')
    [PHUTelem,m] =zipConformingC1_quadratic(PHUTelem,GIFTmesh,patchBoundaries,solIndexCount,p,q);
    
    disp('apply boundary condition')
    [coefSol] = zipConformingC1_BC_circleBiharmonic(PHUTelem,m,p,q);
    
    disp('assign new nodes global and modify c')
    [PHUTelem,sizeBasis,type2Basis] = zipConformingC1_quaddratic_modifyC(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,coefSol,solIndexCount,p,q);
    
    %     disp('localising type 2 basis function')
    %     [PHUTelem] = localizeType2Basis2(PHUTelem,type2Basis);
    %testPlotBasisPhys_GIFTmesh3( PHUTelem,GIFTmesh,p,q,3)
    %testPlotBasisPhys_GIFTmesh(PHUTelem,GIFTmesh,p,q,type2Basis)
    
    %testPlotBasisPhys_GIFTmesh(PHUTelem,GIFTmesh,p,q,type2Basis)
    
    %     for i=1:sizeBasis
    %
    %         testPlotBasisPhys_GIFTmesh3(PHUTelem,GIFTmesh,p,q,i)
    %
    %         pause
    %         close all
    %     end
    
    %%%%%% analysis part ===========================================
    % Material properties
    E  = 10000000;
    nu = 0.3;
    t  = 0.02; % thickness
    
    % Boundary condition
    q0  = -1.;  % distributed force
    clamped = 1;
    
    % constitutive matrix
    D  = E*t^3/(12*(1-nu^2));
    C  = D*[1 nu 0;nu 1 0; 0 0 0.5*(1-nu)];
    
    disp('assembleGalerkin')
    [stiff, rhs] = assembleGalerkinPlate_c1( PHUTelem, GIFTmesh, sizeBasis, p, q, C, q0);
    
    disp('Imposing boundary conditions...')
    [stiff, rhs, bcdof, bcval] = imposeDirichletNeuPHMPCircular(stiff, rhs, PHUTelem, p, q, numPatches,type2Basis);
    
    %toc
    disp('Solving the linear system...')
    %     condest(stiff)
    sol0 = stiff\rhs;
    minZ =PlotDispSquarePlate_c1( PHUTelem, GIFTmesh,  p, q, sol0);
    %     minZ = PlotDispPlate_c1( PHUTelem, GIFTmesh,  p, q, sol0)
    %wbar = minZ*D*1000/q0/W^4
    r=1;
    wbar = minZ
    wext = q0*r^4/64/D
    %pause
    
    %     %%%%%%%%%%%%%%%%%%%%% refinement %%%%%%%%%%%%%%%%%%%%
    if stepCount<numSteps
        numElemU=numElemU*2;
        numElemV=numElemV*2;
        for indexPatch=1:numPatches
            if p==2
                [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmesh_quadratic( p,q,numElemU,numElemV);
            else
                [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmeshGen( p,q,numElemU,numElemV);
            end
            numElem=length(PHUTelem{1});
            % quadList{indexPatch} =1:numElem;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end