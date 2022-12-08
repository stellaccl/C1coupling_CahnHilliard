%exact sol : u(x,y)=(x^2)*((1-x)^2)*(y^2)*((1-y)^2);

restoredefaultpath
close all
clear all
addpaths

tic;

numPatches = 2;
W = 1;
L = 1;
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

numSteps =5;
for stepCount = 1:numSteps
    
    [PHUTelem,  sizeBasis ] = zipConformingNew( PHUTelem, dimBasis, patchBoundaries, p, q);
    
    [PHUTelem ,solIndexCount,patchInfo] = assignSolIndex2D( PHUTelem,patchBoundaries,p,q);
    
    figure
    plotPHTMeshMP(PHUTelem, GIFTmesh)
    
    disp('computing matrix m')
    [m] = zipConformingC1_computeM(PHUTelem,GIFTmesh,patchBoundaries,solIndexCount,p,q);
    
    disp('solving for coef sol')
    [coefSol] = zipConformingC1_boundaryCondition_biharmonic2Patch(PHUTelem,m,p,q );
    
    disp('modifying c')
    [ PHUTelem,sizeBasis,type2Basis] = zipConformingC1_modifyC(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,coefSol,p,q);
    
    disp('localising type 2 basis function')
    [PHUTelem] = localizeType2Basis2(PHUTelem,type2Basis);
    
    % testPlotBasisPhys_GIFTmesh( PHUTelem,GIFTmesh,p,q,type2Basis)
    
    disp('assembleGalerkin')
   
    [stiff,rhs] = assembleBiharmonic_c1_3(PHUTelem,GIFTmesh,sizeBasis,p,q);
    disp('Imposing boundary conditions...')
    [ stiff, rhs, bcdof, bcval] = imposeDirichletBiharmonic_c1(stiff, rhs, PHUTelem, p, q ,type2Basis);
    
    disp('Solving the linear system...')
    sol0 = stiff\rhs;
    
   PlotDispPlate_c1(PHUTelem, GIFTmesh,  p, q, sol0);
    %   PlotDispDiffPlate_c1( PHUTelem, GIFTmesh,  p, q, sol0)
    
    size(sol0)
    [l2relerr] =calcErrorNormsBiharmonic_square(PHUTelem, GIFTmesh,  p, q, sol0 )
    
    if stepCount<numSteps
        for patchIndex=1:numPatches
            quadRef{patchIndex} = 1:size(quadList{patchIndex},1);
            [quadList{patchIndex}, PHUTelem{patchIndex}, dimBasis(patchIndex)] = refineMesh(quadRef{patchIndex}, quadList{patchIndex}, PHUTelem{patchIndex}, p, q, dimBasis(patchIndex));
        end
    end
    
end