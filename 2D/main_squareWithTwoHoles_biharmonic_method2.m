close all
clear all
restoredefaultpath
addpaths
tic
numPatches = 13;
L=1;
W=1;
%GIFTmesh = init2DGeometryGIFTMP('squareWithHole',W,L,numPatches);
%GIFTmesh = init2DGeometryGIFTMP('LShaped_bracket',W,L,numPatches);
GIFTmesh = init2DGeometryGIFTMP('squareWithTwoHoles',W,L,numPatches);
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

% figure
% plotPHTMeshMP(PHUTelem, GIFTmesh)
% axis equal
% pause


patchBoundaries ={1,2,1,3;...
    2,3,4,2;...
    3,4,4,2;...
    4,5,3,1;...
    5,6,3,1;...
    6,7,2,4;...
    [1,7],8,[3,2],[1,4];...
    8,9,2,4;...
    9,10,2,4;...
    10,11,1,3;...
    11,12,1,3;...
    [2,12],13,[2,4],[4,2]};

%%%%  unfinished

numSteps =4;
for stepCount = 1:numSteps
    
    [PHUTelem, dimBasis, quadList ] = checkConforming_IGAPack( PHUTelem, dimBasis, patchBoundaries, p, q, quadList );
    [PHUTelem,  sizeBasis ] = zipConformingNew( PHUTelem, dimBasis, patchBoundaries, p, q);
    %
    %     %     figure
    %     %     plotPHTMeshMP(PHUTelem, GIFTmesh)
    %
    disp('assign sol Index')
    [PHUTelem,solIndexCount] = assignSolIndex2D_quadratic(PHUTelem,patchBoundaries,sizeBasis,p,q);
    % %
    %     figure
    %     plotPHTMesh_solIndex(PHUTelem,GIFTmesh,p)
    %     title('sol Index')
    %     pause
    disp('compute matrix m')
    [PHUTelem,m] =zipConformingC1_quadratic(PHUTelem,GIFTmesh,patchBoundaries,solIndexCount,p,q);
    
    disp('apply boundary condition')
    [coefSol] = zipConformingC1_boundaryCondition_squareWith2Holes( PHUTelem,m,p,q );
    
    disp('assign new nodes global and modify c')
    [PHUTelem,sizeBasis,type2Basis] = zipConformingC1_quaddratic_modifyC(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,coefSol,solIndexCount,p,q);
    %
    %
    disp('localising type 2 basis function')
    [PHUTelem] = localizeType2Basis2(PHUTelem,type2Basis);
    
    type2Basis
    %
    %     for i=type2Basis
    %
    %         %         testPlotBasisPhys_GIFTmesh3(PHUTelem,GIFTmesh,p,q,i)
    %         testPlotBasisPhys_GIFTmesh_flat(PHUTelem,GIFTmesh,p,q,i)
    %         pause
    %         close all
    %     end
    
    %
    disp('assembleGalerkin')
    [stiff, rhs ] = assembleBiharmonic_squareWith2Holes(PHUTelem,GIFTmesh,sizeBasis,p,q);
    %
    %     %%% need to change
    disp('Imposing boundary conditions...')
    [stiff, rhs, bcdof, bcval] = imposeDirichletBiharmonic_squareWith2Holes(stiff, rhs, PHUTelem,GIFTmesh, p, q ,type2Basis);
    %
    disp('Solving the linear system...')
    sol0 = stiff\rhs;
    size(sol0)
    
    [l2relerr] =calcErrorNormsBiharmonic_squareWith2Holes(PHUTelem, GIFTmesh,  p, q, sol0 )
    
%     PlotDispPlate_c1( PHUTelem, GIFTmesh,  p, q, sol0)
%     title(['sol step : ',num2str(stepCount),', error : ',num2str(l2relerr),', numType2Basis : ',num2str(length(type2Basis))])
%     view(0,90)
%     axis tight
%     axis equal
%     colorbar
%   
%    
%     PlotSolDiff_squareWith2Holes( PHUTelem, GIFTmesh,  p, q, sol0)
%     title('diff between sol and exact sol')
%     view(0,90)
%     axis tight
%     axis equal
%     colorbar
%     pause
%     
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
    %
end

%PlotExactSol_squareWith2Holes( PHUTelem, GIFTmesh,  p, q, sol0)
%title('exact sol')