close all
clear all
restoredefaultpath
addpaths
tic
numPatches = 13;
L=1;
W=1;

GIFTmesh = init2DGeometryGIFTMP('squareWithTwoHoles',W,L,numPatches);
p=3;
q=3;


numElemU=2;
numElemV=2;

dimBasis = zeros(1, numPatches);
PHUTelem = cell(numPatches, 1);
quadList = cell(numPatches, 1);

dimBasisDE = zeros(1, numPatches);
PHUTelemDE = cell(numPatches, 1);
quadListDE = cell(numPatches, 1);

if p==2
    
    for indexPatch=1:numPatches
        [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmesh_quadratic( p,q,numElemU,numElemV);
        numElem=length(PHUTelem{indexPatch});
        quadList{indexPatch} =1:numElem;
        
    end
    
    for indexPatch=1:numPatches
        [PHUTelemDE{indexPatch}, dimBasisDE(indexPatch)] = initPHTmeshGen( p+1,q+1,numElemU,numElemV);
        numElem=length(PHUTelemDE{indexPatch});
        quadListDE{indexPatch} =1:numElem;
    end
    
else
    
    for indexPatch=1:numPatches
        [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmeshGen( p,q,numElemU,numElemV);
        numElem=length(PHUTelem{indexPatch});
        quadList{indexPatch} =1:numElem;
        
        [PHUTelemDE{indexPatch}, dimBasisDE(indexPatch)] = initPHTmeshGen( p+1,q+1,numElemU,numElemV);
        numElem=length(PHUTelemDE{indexPatch});
        quadListDE{indexPatch} =1:numElem;
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
    
    [PHUTelemDE, dimBasisDE, quadListDE ] = checkConforming_IGAPack( PHUTelemDE, dimBasisDE, patchBoundaries, p+1, q+1, quadListDE );
    [PHUTelemDE,  sizeBasisDE ] = zipConformingNew( PHUTelemDE, dimBasisDE, patchBoundaries, p+1, q+1);

    
    %     %     figure
    %     %     plotPHTMeshMP(PHUTelem, GIFTmesh)
    %
    disp('assign sol Index')
    [PHUTelem,~] = assignSolIndex2D_DE_2(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,p,q);
    [PHUTelemDE,solIndexCount] = assignSolIndex2D_DE_2(PHUTelemDE,GIFTmesh,patchBoundaries,sizeBasisDE,p+1,q+1);

    disp('remove type2 basis function')
    [PHUTelem,sizeBasis] = removeType2Basis_degElev(PHUTelem,p,q ,sizeBasis);
 
    disp('compute matrix m')
    [PHUTelemDE,m] =zipConformingC1_quadratic(PHUTelemDE,GIFTmesh,patchBoundaries,solIndexCount,p+1,q+1);
  
    disp('apply boundary condition')
    [coefSol] = zipConformingC1_boundaryCondition_squareWith2Holes( PHUTelemDE,m,p+1,q+1 );
    
    % %
    %     figure
    %     plotPHTMesh_solIndex(PHUTelem,GIFTmesh,p)
    %     title('sol Index')
    %     pause
   
    disp('assign type 2 basis')
    [PHUTelem,sizeBasis,type2Basis] = assignC1basis_degElev(PHUTelem,PHUTelemDE,sizeBasis,coefSol,solIndexCount,p,q);
    
    disp('assign additional basis')
    [PHUTelem,PHUTelemDE,sizeBasis ] = assignAdditionalNodes_degElev( PHUTelem,PHUTelemDE,sizeBasis,p,q);
  
    
    %     for i=type2Basis
    %
    %         %         testPlotBasisPhys_GIFTmesh3(PHUTelem,GIFTmesh,p,q,i)
    %         testPlotBasisPhys_GIFTmesh_flat(PHUTelem,GIFTmesh,p,q,i)
    %         pause
    %         close all
    %     end
    
    %
    disp('assembleGalerkin')
    [stiff, rhs ] = assembleBiharmonic_squareWith2Holes_degElev(PHUTelem,GIFTmesh,sizeBasis,p,q);
    %
    %     %%% need to change
    disp('Imposing boundary conditions...')
    [stiff, rhs, bcdof, bcval] = imposeDirichletBiharmonic_squareWith2Holes_degElev(stiff, rhs, PHUTelem,PHUTelemDE,GIFTmesh, p, q ,type2Basis);
    %
    disp('Solving the linear system...')
    sol0 = stiff\rhs;
    size(sol0)
    
    [l2relerr] =calcErrorNormsBiharmonic_squareWith2Holes_degElev(PHUTelem, GIFTmesh,  p, q, sol0 );
    
    
     PlotDispPlate_c1_degElev(PHUTelem, GIFTmesh,  p, q, sol0);
    title(['sol step : ',num2str(stepCount),', error : ',num2str(l2relerr),', numType2Basis : ',num2str(length(type2Basis))])
    view(0,90)
    axis tight
    axis equal
    colorbar
    drawnow
   
%     PlotSolDiff_squareWith2Holes( PHUTelem, GIFTmesh,  p, q, sol0)
%     title('diff between sol and exact sol')
%     view(0,90)
%     axis tight
%     axis equal
%     colorbar
%     pause
    
    %%%%%%%%%%%%%%%%%%%%% refinement %%%%%%%%%%%%%%%%%%%%
    if stepCount<numSteps
        numElemU=numElemU*2;
        numElemV=numElemV*2;
        for indexPatch=1:numPatches
            if p==2
                [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmesh_quadratic( p,q,numElemU,numElemV);
                [PHUTelemDE{indexPatch}, dimBasisDE(indexPatch)] = initPHTmeshGen( p+1,q+1,numElemU,numElemV);
            else
                [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmeshGen( p,q,numElemU,numElemV);
                [PHUTelemDE{indexPatch}, dimBasisDE(indexPatch)] = initPHTmeshGen( p+1,q+1,numElemU,numElemV);
            end
            numElem=length(PHUTelem{1});
            quadList{indexPatch} =1:numElem;
            numElem=length(PHUTelemDE{1});
            quadListDE{indexPatch} =1:numElem;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
end

%PlotExactSol_squareWith2Holes( PHUTelem, GIFTmesh,  p, q, sol0)
%title('exact sol')