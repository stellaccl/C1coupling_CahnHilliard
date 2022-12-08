close all
clear all
restoredefaultpath
addpaths
tic
numPatches = 4;
L=1;
W=1;

GIFTmesh = init2DGeometryGIFTMP('torus2',W,L,numPatches);

p=3;
q=3;


numElemU=2;
numElemV=2;

PHUTelem = cell(numPatches, 1);
dimBasis = zeros(1, numPatches);
quadList = cell(numPatches, 1);
if p==2
    for indexPatch=1:numPatches
        [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmesh_quadratic( p,q,numElemU,numElemV);
        numElem=length(PHUTelem{indexPatch});
        quadList{indexPatch} =1:numElem;
    end
else
    for indexPatch=1:numPatches
        [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmeshGen( p,q,numElemU,numElemV);
        numElem=length(PHUTelem{indexPatch});
        quadList{indexPatch} =1:numElem;
    end
end



patchBoundaries ={1,2,1,3;...
    2,3,1,3;...
    [1,3],4,[3,1],[1,3]};


numSteps =3;
for stepCount = 1:numSteps
    
    [PHUTelem, dimBasis, quadList ] = checkConforming_IGAPack( PHUTelem, dimBasis, patchBoundaries, p, q, quadList );
    [PHUTelem,  sizeBasis ] = zipConformingNew( PHUTelem, dimBasis, patchBoundaries, p, q);
    
    %     figure
    %     plotPHTMeshMP(PHUTelem, GIFTmesh)
    
    
    disp('assign sol Index')
    [PHUTelem,solIndexCount] = assignSolIndex2D_quadratic(PHUTelem,patchBoundaries,sizeBasis,p,q);
    
%     
%     figure
%     plotPHTMesh_solIndex(PHUTelem,GIFTmesh,p)
%     title('sol Index')
    
    disp('compute matrix m')
    [PHUTelem,m] =zipConformingC1_quadratic(PHUTelem,GIFTmesh,patchBoundaries,solIndexCount,p,q);
    
    
    disp('apply boundary condition')
    [coefSol] = zipConformingC1_quadratic_boundaryCondition_annulus(PHUTelem,m,p,q);
    
    disp('assign new nodes global and modify c')
    [PHUTelem,sizeBasis,type2Basis] = zipConformingC1_quaddratic_modifyC(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,coefSol,solIndexCount,p,q);
    
    
    disp('localising type 2 basis function')
    [PHUTelem] = localizeType2Basis2(PHUTelem,type2Basis);
    
    disp('assembleGalerkin')
    % [stiff,rhs] = assembleBiharmonicAnnulus_c1(PHUTelem,GIFTmesh,sizeBasis,p,q);
    
    
    [stiff, rhs ] = assembleBiharmonicAnnulus_2(PHUTelem,GIFTmesh,sizeBasis,p,q);
    
    
    disp('Imposing boundary conditions...')
    [stiff, rhs, bcdof, bcval] = imposeDirichletBiharmonic_annulus(stiff, rhs, PHUTelem,GIFTmesh, p, q ,type2Basis);
    
    disp('Solving the linear system...')
    sol0 = stiff\rhs;
    size(sol0)
    [l2relerr] =calcErrorNormsBiharmonic_annulus(PHUTelem, GIFTmesh,  p, q, sol0 )
    
    PlotDispPlate_c1( PHUTelem, GIFTmesh,  p, q, sol0)
    axis equal
    
    PlotDispPlate_annulus_diff( PHUTelem, GIFTmesh,  p, q, sol0)
%    pause
    
    %%%%%%%%%%%%%%%%%%%%% refinement %%%%%%%%%%%%%%%%%%%%
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
            quadList{indexPatch} =1:numElem;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

PlotDispPlate_c1( PHUTelem, GIFTmesh,  p, q, sol0)
axis equal