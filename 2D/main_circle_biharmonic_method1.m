%5patch circle example
close all
clear all
restoredefaultpath
addpaths
tic
numPatches = 5;
L=1;
W=1;

GIFTmesh = init2DGeometryGIFTMP('5PatchCircle',W,L,numPatches);

dimBasis = zeros(1, numPatches);
p=4;
q=4;

numberElementsU=2;
numberElementsV=2;

PHUTelem = cell(numPatches, 1);
for i=1:numPatches
    [PHUTelem{i}, dimBasis(i)] = initPHTmesh( p,q,numberElementsU,numberElementsV );
    quadList{i} = 2:5;
end

figure
plotPHTMeshMP(PHUTelem, GIFTmesh)

patchBoundaries = {1,2,2,2;...
    2,3,4,2;...
    [1,3],4,[4,4],[2,4];...
    [1,2,3,4],5,[3,3,3,3],[2,3,4,1]};

[PHUTelem,  sizeBasis ] = zipConformingNew( PHUTelem, dimBasis, patchBoundaries, p, q);

figure
plotPHTMeshMP( PHUTelem, GIFTmesh )

disp('assign sol Index')
[PHUTelem,solIndexCount] = assignSolIndex2D(PHUTelem,patchBoundaries,p,q);

figure
plotPHTMesh_nodesGlobal(PHUTelem, GIFTmesh,p)
title('nodes Global')

figure
plotPHTMesh_solIndex( PHUTelem,GIFTmesh,p)
title('sol Index')

disp('computing matrix m')
[m] = zipConformingC1_computeM(PHUTelem,GIFTmesh,patchBoundaries,solIndexCount,p,q);

disp('solving for coef sol')
[coefSol] = zipConformingC1_boundaryCondition_circleBiharmonic(PHUTelem,m,p,q);

disp('modifying c')
[PHUTelem,sizeBasis,type2Basis] = zipConformingC1_modifyC(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,coefSol,p,q);

[PHUTelem] = localizeType2Basis2(PHUTelem,type2Basis);

%testPlotBasisPhys_GIFTmesh( PHUTelem,GIFTmesh,p,q,sizeBasis,numType2Basis)
disp('assembleGalerkin')
%[stiff,rhs] = assembleBiharmonicCircle_c1(PHUTelem,GIFTmesh,sizeBasis,p,q);
[stiff,rhs] = assembleBiharmonicCircle_2(PHUTelem,GIFTmesh,sizeBasis,p,q);
disp('Imposing boundary conditions...')
[stiff, rhs, bcdof, bcval] = imposeDirichletBiharmonic_circle(stiff, rhs, PHUTelem,GIFTmesh, p, q ,type2Basis);

disp('Solving the linear system...')
sol0 = stiff\rhs;
%     [l2relerr] =calcErrorNormsBiharmonic_annulus(PHUTelem, GIFTmesh,  p, q, sol0 )


PlotDispPlate_c1( PHUTelem, GIFTmesh,  p, q, sol0)
axis equal

