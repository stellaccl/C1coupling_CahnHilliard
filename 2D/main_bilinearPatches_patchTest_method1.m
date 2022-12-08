%2patch generic & non generic example
close all
clear all

restoredefaultpath
addpaths

numPatches = 2;
GIFTmesh = init2DGeometryGIFTMP('bilinear_rectangle');

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

patchBoundaries = {1, 2, 2, 4};

%define the number of refinement steps
numSteps = 2;
for stepCount = 1:numSteps
    
    [PHUTelem, dimBasis, quadList ] = checkConforming_IGAPack( PHUTelem, dimBasis, patchBoundaries, p, q, quadList );
    [PHUTelem, sizeBasis ] = zipConformingNew( PHUTelem, dimBasis, patchBoundaries, p, q);
    
    %     figure
    %     plotPHTMeshMP( PHUTelem, GIFTmesh )
    %     pause
    %     figure
    %     plotPHTMesh_nodesGlobal(PHUTelem, GIFTmesh,p)
    %     title('nodes Global')
    
    disp('assign sol Index')
    [PHUTelem,solIndexCount] = assignSolIndex2D(PHUTelem,patchBoundaries,p,q);
    
    %     figure
    %     plotPHTMesh_solIndex( PHUTelem,GIFTmesh,p)
    %     title('sol Index')
    
    
    disp('computing matrix m')
    [m] = zipConformingC1_computeM(PHUTelem,GIFTmesh,patchBoundaries,solIndexCount,p,q);
    
    disp('solving for coef sol')
    [coefSol] = zipConformingC1_boundaryCondition_bilinearPatch_patchTest(PHUTelem,m,p,q);
    
    disp('modifying c')
    [PHUTelem,sizeBasis,type2Basis] = zipConformingC1_modifyC(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,coefSol,p,q);
    
    [PHUTelem] = localizeType2Basis2(PHUTelem,type2Basis);
    
    Emod=1e5;
    nu=0;
    Cmat = Emod/((1+nu)*(1-2*nu))*[1-nu, nu, 0; nu, 1-nu, 0; 0, 0, (1-2*nu)/2];
    bound_disp = 0.1;
    
    disp('Assembling the linear system...')
    [ stiff, rhs ] = assembleGalerkinSysGIFTMP_modifiedC( PHUTelem, GIFTmesh, sizeBasis, p, q, Cmat);
    
    disp('Imposing boundary conditions...')
    [ stiff, rhs, bcdof, bcval] = imposeDirichletPTestPHUT_billinear2patch(stiff, rhs, PHUTelem, p, q, bound_disp,patchBoundaries);
    condest(stiff)
    
    disp('Solving the linear system...')
    sol0 = stiff\rhs;
    
    disp('Plotting the solution...')
    plotSolPHUTElasticVM2_GIFTmesh( sol0, GIFTmesh,PHUTelem, p, q, Cmat)
    
    %uniform refinement
    for patchIndex=1:numPatches
        quadRef{patchIndex} = 1:size(quadList{patchIndex},1);
        [quadList{patchIndex}, PHUTelem{patchIndex}, dimBasis(patchIndex)] = refineMesh(quadRef{patchIndex}, quadList{patchIndex}, PHUTelem{patchIndex}, p, q, dimBasis(patchIndex));
    end
%    pause
end
