% for all p

close all
clear all
restoredefaultpath

addpaths
numPatches = 3;
numElemU=6;
numElemV=6;
W = 100;
L = 100;
p=4;
q=4;

GIFTmesh = init2DGeometryGIFTMP('threePatch');

dimBasis = zeros(1, numPatches);
PHUTelem = cell(numPatches, 1);

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



patchBoundaries = {1,2,3,1;...
    [1,2],3,[4,4],[1,2]};

numSteps =3;
for stepCount = 1:numSteps
    
    figure
    plotPHTMeshMP(PHUTelem, GIFTmesh)
    
    [PHUTelem,  sizeBasis ] = zipConformingNew( PHUTelem, dimBasis, patchBoundaries, p, q);
    
    disp('assign sol Index')
    [PHUTelem,solIndexCount] = assignSolIndex2D_quadratic(PHUTelem,patchBoundaries,sizeBasis,p,q);
    
    % figure
    % plotPHTMesh_nodesGlobal( PHUTelem,GIFTmesh,p)
    % title('nodes Global')
    % pause
    %
    % figure
    % plotPHTMesh_solIndex(PHUTelem,GIFTmesh,p)
    % title('sol Index')
    % pause
    
    % figure
    % plotPHTMesh_solIndex( PHUTelem,GIFTmesh,p)
    % title('sol Index after assign neighbor')
    % pause
    
    disp('compute matrix m')
    [PHUTelem,m] =zipConformingC1_quadratic(PHUTelem,GIFTmesh,patchBoundaries,solIndexCount,p,q);
   
    disp('apply boundary condition')
    [coefSol] = zipConformingC1_p2_BC_3PatchBiharmonic(PHUTelem,m,p,q);
   
    disp('assign new nodes global and modify c')
    [PHUTelem,sizeBasis,type2Basis] = zipConformingC1_quaddratic_modifyC(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,coefSol,solIndexCount,p,q);
    
%     type2Basis
%     pause
    disp('localising type 2 basis function')
    [PHUTelem] = localizeType2Basis2(PHUTelem,type2Basis);
    
%     testPlotBasisPhys_GIFTmesh(PHUTelem,GIFTmesh,p,q,type2Basis)
    
    disp('assembleGalerkin')
    [stiff,rhs] = assembleBiharmonic_c1_3(PHUTelem,GIFTmesh,sizeBasis,p,q);
    
    disp('Imposing boundary conditions...')
    [ stiff, rhs, bcdof, bcval ] = imposeDirichlet_3patch_biharmonic(stiff, rhs, PHUTelem, p, q, numPatches,type2Basis);
    
    disp('Solving the linear system...')
    sol0 = stiff\rhs;
    
    PlotDispPlate_c1(PHUTelem, GIFTmesh,  p, q, sol0);
    axis equal
    
    size(sol0)
    [l2relerr] =calcErrorNormsBiharmonic_square(PHUTelem, GIFTmesh, p, q, sol0 )
    
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