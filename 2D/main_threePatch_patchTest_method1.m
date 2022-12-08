close all
clear all
restoredefaultpath

addpaths
numPatches = 3;
numberElementsU=4;
numberElementsV=4;
p=3;
q=3;

GIFTmesh = init2DGeometryGIFTMP('threePatch');

dimBasis = zeros(1, numPatches);
PHUTelem = cell(numPatches, 1);
for i=1:numPatches
    [PHUTelem{i}, dimBasis(i)] = initPHTmeshGen(p,q,numberElementsU,numberElementsV);
%    quadList{indexPatch} = 1:numberElementsU*numberElementsV;
end

figure
plotPHTMeshMP(PHUTelem, GIFTmesh)
%pause

patchBoundaries = {1,2,3,1;...
    [1,2],3,[4,4],[1,2]};


[PHUTelem,  sizeBasis ] = zipConformingNew( PHUTelem, dimBasis, patchBoundaries, p, q);

disp('assigning solution index')
[PHUTelem,solIndexCount] = assignSolIndex2D(PHUTelem,patchBoundaries,p,q);

disp('computing matrix m')
[m] = zipConformingC1_computeM(PHUTelem,GIFTmesh,patchBoundaries,solIndexCount,p,q);

disp('solving for coef sol')
[coefSol] = zipConformingC1_boundaryCondition_threePatch(PHUTelem,m,p,q );

disp('modifying c')
[ PHUTelem,sizeBasis,type2Basis] = zipConformingC1_modifyC(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,coefSol,p,q);

% figure
% plotPHTMesh_nodesGlobal2( PHUTelem,GIFTmesh,type2Basis,p);
% title('nodes Global 2')
% pause


disp('localising type 2 basis function')
[PHUTelem] = localizeType2Basis2(PHUTelem,type2Basis);

%testPlotBasisPhys_GIFTmesh(PHUTelem,GIFTmesh,p,q,sizeBasis,numType2Basis)

%analysis part
Emod=1e5;
nu=0;
P = 1e4;
Cmat = Emod/((1+nu)*(1-2*nu))*[1-nu, nu, 0; nu, 1-nu, 0; 0, 0, (1-2*nu)/2];
bound_disp = 0.1;


disp('Assembling the linear system...')
[ stiff, rhs ] = assembleGalerkinSysGIFTMPC1( PHUTelem, GIFTmesh, sizeBasis, p, q, Cmat );
condest(stiff)


disp('Imposing boundary conditions...')
[ stiff, rhs, bcdof, bcval] = imposeDirichletNeu_threePatch(stiff, rhs, PHUTelem,GIFTmesh, p, q,P,type2Basis);


disp('Solving the linear system...')
sol0 = stiff\rhs;

disp('Plotting the solution...')
plotSolPHUTElasticVM2_GIFTmesh(sol0, GIFTmesh,PHUTelem, p, q, Cmat)

