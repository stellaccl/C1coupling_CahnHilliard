close all
clear all
restoredefaultpath
addpaths
tic
numPatches = 4;
L=1;
W=1;
%GIFTmesh = init2DGeometryGIFTMP('squareWithHole',W,L,numPatches);
GIFTmesh = init2DGeometryGIFTMP('LShaped_bracket',W,L,numPatches);
p=4;
q=4;


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


% 
% figure
%    plotPHTMeshMP(PHUTelem, GIFTmesh)
%    axis equal
% pause

patchBoundaries ={1,2,3,1;...
    2,3,3,1;...
    [1,3],4,[1,3],[3,1]};

numSteps =3;
for stepCount = 1:numSteps
    
    [PHUTelem, dimBasis, quadList ] = checkConforming_IGAPack( PHUTelem, dimBasis, patchBoundaries, p, q, quadList );
    [PHUTelem,  sizeBasis ] = zipConformingNew( PHUTelem, dimBasis, patchBoundaries, p, q);
    
    %     figure
    %     plotPHTMeshMP(PHUTelem, GIFTmesh)
    
    disp('assign sol Index')
    [PHUTelem,solIndexCount] = assignSolIndex2D_quadratic(PHUTelem,patchBoundaries,sizeBasis,p,q);
   
%     figure
%     plotPHTMesh_solIndex(PHUTelem,GIFTmesh,p)
%     title('sol Index')
%    
    disp('compute matrix m')
    [PHUTelem,m] =zipConformingC1_quadratic(PHUTelem,GIFTmesh,patchBoundaries,solIndexCount,p,q);

    disp('apply boundary condition')
    [coefSol] = zipConformingC1_boundaryCondition_squareWithHole( PHUTelem,m,p,q );
    
    disp('assign new nodes global and modify c')
    [PHUTelem,sizeBasis,type2Basis] = zipConformingC1_quaddratic_modifyC(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,coefSol,solIndexCount,p,q);
    
    
    disp('localising type 2 basis function')
    [PHUTelem] = localizeType2Basis2(PHUTelem,type2Basis);
    
    
%     for i=type2Basis
%         
% %         testPlotBasisPhys_GIFTmesh3(PHUTelem,GIFTmesh,p,q,i)
%          testPlotBasisPhys_GIFTmesh_flat(PHUTelem,GIFTmesh,p,q,i)
%         pause
%         close all
%     end
    
    
   
    disp('assembleGalerkin')
    [stiff, rhs ] = assembleBiharmonic_squareWithHole(PHUTelem,GIFTmesh,sizeBasis,p,q);
    
    %%% need to change
    disp('Imposing boundary conditions...')
    [stiff, rhs, bcdof, bcval] = imposeDirichletBiharmonic_squareWithHole(stiff, rhs, PHUTelem,GIFTmesh, p, q ,type2Basis);
    
    disp('Solving the linear system...')
    sol0 = stiff\rhs;
    size(sol0)
    [l2relerr] =calcErrorNormsBiharmonic_squareWithHole(PHUTelem, GIFTmesh,  p, q, sol0 )
    
    PlotDispPlate_c1( PHUTelem, GIFTmesh,  p, q, sol0)
    view(0,90)
    axis tight
    axis equal
    colorbar
    pause
    
    %%% need to change
%     PlotDispPlate_annulus_diff( PHUTelem, GIFTmesh,  p, q, sol0)
%     pause
    
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

