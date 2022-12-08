restoredefaultpath
close all
clear all
addpaths

tic;

numPatches = 2;
W = 1;
L = 1;
%GIFTmesh = init2DGeometryGIFTMP('rectangle', L, W, numPatches);
GIFTmesh = init2DGeometryGIFTMP('2patch_curvedCommonBoundary', L, W, numPatches);
p=3;
q=3;

numElemU=3;
numElemV=3;

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

patchBoundaries = {1,2,2,4};

numSteps =4;
for stepCount = 1:numSteps
    
    figure
    plotPHTMeshMP(PHUTelem, GIFTmesh)
    %pause
    
    [PHUTelem,  sizeBasis ] = zipConformingNew( PHUTelem, dimBasis, patchBoundaries, p, q);
    
    disp('assign sol Index')
    [PHUTelem,solIndexCount] = assignSolIndex2D_quadratic(PHUTelem,patchBoundaries,sizeBasis,p,q);
    
    disp('compute matrix m')
    [PHUTelem,m] =zipConformingC1_quadratic(PHUTelem,GIFTmesh,patchBoundaries,solIndexCount,p,q);
   % [PHUTelem,m] =zipConformingC1_method3(PHUTelem,GIFTmesh,patchBoundaries,solIndexCount,p,q);
    disp('apply boundary condition')
    [coefSol] = zipConformingC1_p2_BC_2patchSquareBiharmonic(PHUTelem,m,p,q);
    
    disp('assign new nodes global and modify c')
    [PHUTelem,sizeBasis,type2Basis] = zipConformingC1_quaddratic_modifyC(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,coefSol,solIndexCount,p,q);
    
    %     figure
    %     plotPHTMesh_solIndex( PHUTelem,GIFTmesh,p)
    %     title('sol Index')
    %     pause
    %
    %     figure
    %     plotPHTMesh_nodesGlobal( PHUTelem,GIFTmesh,p)
    %     title('nodes Global')
    %     pause
    
    disp('localising type 2 basis function')
    [PHUTelem] = localizeType2Basis2(PHUTelem,type2Basis);
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
    
    Emod=1e5;
    nu=0;
    Cmat = Emod/((1+nu)*(1-2*nu))*[1-nu, nu, 0; nu, 1-nu, 0; 0, 0, (1-2*nu)/2];
    bound_disp = 0.1;
    
    disp('Assembling the linear system...')
    [stiff,rhs] = assembleBiharmonic_c1_3(PHUTelem,GIFTmesh,sizeBasis,p,q);
    disp('Imposing boundary conditions...')
    [ stiff, rhs, bcdof, bcval] = imposeDirichletBiharmonic_c1(stiff, rhs, PHUTelem, p, q ,type2Basis);
    
    disp('Solving the linear system...')
    sol0 = stiff\rhs;
    
    figure
    PlotDispPlate_c1(PHUTelem, GIFTmesh,  p, q, sol0);
    axis equal
    
    size(sol0)
    [l2relerr] =calcErrorNormsBiharmonic_square(PHUTelem, GIFTmesh,  p, q, sol0)
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

PlotDispPlate_c1(PHUTelem, GIFTmesh,  p, q, sol0);
axis equal