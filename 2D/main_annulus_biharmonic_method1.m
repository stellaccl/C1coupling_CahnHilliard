close all
clear all
restoredefaultpath
addpaths
tic
numPatches = 4;
L=1;
W=1;

GIFTmesh = init2DGeometryGIFTMP('torus2',W,L,numPatches);

dimBasis = zeros(1, numPatches);
p=3;
q=3;

numElemU=1;
numElemV=1;
PHUTelem = cell(numPatches, 1);
for i=1:numPatches
    [PHUTelem{i}, dimBasis(i)] = initPHTmesh(p,q,numElemU,numElemV);
    quadList{i} = 2:5;
end


patchBoundaries ={1,2,1,3;...
    2,3,1,3;...
    [1,3],4,[3,1],[1,3]};


numSteps =2;
for stepCount = 1:numSteps
    
%     figure
%     plotPHTMeshMP(PHUTelem, GIFTmesh)
    
    disp('zipConforming...')
    [PHUTelem,sizeBasis] = zipConformingNew( PHUTelem, dimBasis, patchBoundaries, p, q);
    
    disp('assign sol Index')
    [PHUTelem,solIndexCount] = assignSolIndex2D(PHUTelem,patchBoundaries,p,q);
    
    disp('computing matrix m')
    [m] = zipConformingC1_computeM(PHUTelem,GIFTmesh,patchBoundaries,solIndexCount,p,q);
    
    disp('solving for coef sol')
    [coefSol] = zipConformingC1_boundaryCondition_annulusBiharmonic(PHUTelem,m,p,q);
    
    disp('modifying c')
    [PHUTelem,sizeBasis,type2Basis] = zipConformingC1_modifyC(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,coefSol,p,q);
    
    disp('localising type 2 basis function')
    [PHUTelem] = localizeType2Basis2(PHUTelem,type2Basis);
    
    disp('assembleGalerkin')
    [stiff,rhs] = assembleBiharmonicAnnulus_c1(PHUTelem,GIFTmesh,sizeBasis,p,q);
    
    disp('Imposing boundary conditions...')
    [stiff, rhs, bcdof, bcval] = imposeDirichletBiharmonic_annulus(stiff, rhs, PHUTelem,GIFTmesh, p, q ,type2Basis);
  
    disp('Solving the linear system...')
    sol0 = stiff\rhs;
    [l2relerr] =calcErrorNormsBiharmonic_annulus(PHUTelem, GIFTmesh,  p, q, sol0 )

     
    PlotDispPlate_c1( PHUTelem, GIFTmesh,  p, q, sol0)
   axis equal

    if stepCount<numSteps
        for patchIndex=1:numPatches
            quadRef{patchIndex} = 1:size(quadList{patchIndex},1);
            [quadList{patchIndex}, PHUTelem{patchIndex}, dimBasis(patchIndex)] = refineMesh(quadRef{patchIndex}, quadList{patchIndex}, PHUTelem{patchIndex}, p, q, dimBasis(patchIndex));
        end
    end
    
end

PlotExactSol_annulusBiharmonic(PHUTelem, GIFTmesh,  p, q)