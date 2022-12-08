% PlateHole example
%plotErrorPHTElasticPH_modifiedC
%plotSolPHTElasticMP_modifiedC
%Use ComputeMatixM

close all
clear all

restoredefaultpath
addpaths

numPatches = 2;
p=4;
q=4;

GIFTmesh = init2DGeometryGIFTMP('plate_hole');
numElemU=1;
numElemV=1;
dimBasis = zeros(1, numPatches);
PHUTelem = cell(numPatches, 1);
for indexPatch=1:numPatches
    [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmesh(p,q,numElemU,numElemV);
    quadList{indexPatch} = 2:5;
end
% 
% figure
% plotPHTMeshMP(PHUTelem, GIFTmesh)


patchBoundaries = {1, 2, 2, 4};

%define the number of refinement steps
numSteps =5;
for stepCount = 1:numSteps
    
    
    [PHUTelem,  sizeBasis ] = zipConformingNew( PHUTelem, dimBasis, patchBoundaries, p, q);
    
    figure
    plotPHTMeshMP( PHUTelem, GIFTmesh )
    
    disp('assign sol Index')
    [PHUTelem,solIndexCount] = assignSolIndex2D(PHUTelem,patchBoundaries,p,q);   
   
    disp('zip conforming C1')
    [PHUTelem,sizeBasis,numType2Basis,type2Basis] =zipConformingC1(PHUTelem,GIFTmesh,sizeBasis,patchBoundaries,solIndexCount,p,q);
    %   testPlotBasisPhys_GIFTmesh( PHUTelem,GIFTmesh,p,q,sizeBasis,numType2Basis)
    
    
    %%%======================= analysis part ==============================
    %Material properties
    Emod = 1e5;
    nu = 0;
    
    %elasticity matrix (plane stress)
    Cmat=Emod/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
    
    %the traction at infinity
    tx = 10;
    
    %define radius of the hole and side length
    rad = 1;
    L = 4;
    
    disp('Assembling the linear system...')
    [ stiff, rhs ] = assembleGalerkinSysGIFTMP_modifiedC( PHUTelem, GIFTmesh, sizeBasis, p, q, Cmat);
    
    disp('Imposing boundary conditions...')
    [ stiff, rhs, bcdof, bcval ] = imposeDirichletNeu_plateWithHole(stiff, rhs, PHUTelem, GIFTmesh, p, q, rad, tx);
    
    %toc
    disp('Solving the linear system...')
    stepCount
    
    %     cond(full(stiff))
    %     sol0 = stiff\rhs;
    M1=spdiags(1./sqrt(diag(stiff)),0,size(stiff,1),size(stiff,2));
    LHS = M1*stiff*M1;
    sol0 = LHS\(M1*rhs);
    sol0 = M1*sol0;
    size(sol0)
    %cond(full(LHS))
    
    disp('Computing the errors...')
    [l2relerr, h1relerr]=calcErrorNormsPlateHole_modifiedC( sol0, PHUTelem, GIFTmesh, p, q, Cmat, Emod, nu, rad, tx );
    %     l2relerr
    %     h1relerr
    
    disp('Plotting the errors...')
    plotErrorPHTElasticPlateHoleC1( sol0, PHUTelem, GIFTmesh, p, q, Cmat, rad, tx )
    
    %uniform refinement
    for patchIndex=1:numPatches
        quadRef{patchIndex} = 1:size(quadList{patchIndex},1);
        [quadList{patchIndex}, PHUTelem{patchIndex}, dimBasis(patchIndex)] = refineMesh(quadRef{patchIndex}, quadList{patchIndex}, PHUTelem{patchIndex}, p, q, dimBasis(patchIndex));
    end
end
