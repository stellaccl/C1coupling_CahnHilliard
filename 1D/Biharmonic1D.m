restoredefaultpath
close all
clear
addpath('./utils/')
addpath('./nurbs_toolbox')

tic;

numPatches = 2;
W = 1;
L = 1;

GIFTmesh = init1DGeometryGIFTMP('straightLine');
%
p=2;
numElemU=3;

dimBasis = zeros(1, numPatches);
PHUTelem = cell(numPatches, 1);

if p==2
    for indexPatch=1:numPatches
        [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmesh_quadratic1D( p,numElemU);
    end
else
    
    for indexPatch=1:numPatches
        [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmeshGen_1D( p,numElemU);
    end
end


numSteps =5;
for stepCount = 1:numSteps
    
%     figure
%     plotPHTMeshMP1D(PHUTelem, GIFTmesh)
    %pause
    %zipConformingC0
    [PHUTelem, sizeBasis] =zipConformingC0_1D(PHUTelem);
    
    %figure
%     [colors ] = colorGen([0, 0.5, 0],[1 1 0], sizeBasis);
%     for index=1:sizeBasis
%         
%         plotBasis1D_GIFTmesh_2(PHUTelem,GIFTmesh,p,index,colors(index,:))
%     end
    
    [PHUTelem,solIndexCount] = assignSolIndex1D(PHUTelem,p);
    
    [PHUTelem,type2Basis,sizeBasis] = zipConformingC1_1D(PHUTelem,GIFTmesh,p,sizeBasis,solIndexCount);
    
    %plot C1 basis
%     figure
%     plotPHTMeshMP1D(PHUTelem, GIFTmesh)
%     
%     [colors ] = colorGen([0, 0.5, 0],[1 1 0], sizeBasis);
%     
%     for index=1:sizeBasis
%         % plotBasis1D_modifiedC(PHUTelem,GIFTmesh,p,index)
%         plotBasis1D_modifiedC_paper(PHUTelem,GIFTmesh,p,index,colors(index,:))
%     end
   
    
    [stiff, rhs ] = assembleBiharmonic_c1_1D(PHUTelem,GIFTmesh, sizeBasis, p );
    [stiff, rhs, bcdof, bcval] = imposeDirichletBiharmonic_1D(stiff, rhs, PHUTelem, p,type2Basis);
    
    disp('Solving the linear system...')
    sol0 = stiff\rhs;
    size(sol0)
    
    PlotDisp1D(PHUTelem, GIFTmesh, p, sol0)
    plotExactSolBiharmonic1D(PHUTelem, GIFTmesh)
    format long
    [l2relerr] =calErrorNormBiharmonic_1D( PHUTelem, GIFTmesh, p, sol0)
    title(['step',num2str(stepCount),' l2relerr=',num2str(l2relerr)])
    pause
    %%%%%%%%%%%%%%%%%%%%% refinement %%%%%%%%%%%%%%%%%%%%
    if stepCount<numSteps
        numElemU=numElemU*2;
        for indexPatch=1:numPatches
            if p==2
                [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmesh_quadratic1D( p,numElemU);
            else
                [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmeshGen_1D( p,numElemU);
            end
            numElem=length(PHUTelem{1});
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end