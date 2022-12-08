%PlateHole example with degree elevation
%plotErrorPHTElasticPH_modifiedC
%plotSolPHTElasticMP_modifiedC
%Use ComputeMatixM
%refine elementA and B

close all
clear all
restoredefaultpath
addpaths
tic
numPatches = 2;
GIFTmesh = init2DGeometryGIFTMP('plate_hole');
%GIFTmesh = init2DGeometryGIFTMP('bilinear_wedge');
%targetScale = 1 --> uniform refinement
%tagetScale < 1 --> more graded refinement
targetScale = 0.50;
target_rel_error = 1e-5;

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

%initialize the degree elevated mesh
PHUTelemDE = cell(numPatches, 1);
for i=1:numPatches
    [PHUTelemDE{i}, dimBasisDE(i)] = initPHTmesh(p+1,q+1,numElemU,numElemV);
    quadList{i} = 2:5;
end

% figure
% plotPHTMeshMP(PHUTelem, GIFTmesh)

patchBoundaries = {1, 2, 2, 4};

%define the number of refinement steps
numSteps =3;
for stepCount = 1:numSteps
    
    [PHUTelem,  sizeBasis ] = zipConformingNew( PHUTelem, dimBasis, patchBoundaries, p, q);
    [PHUTelemDE,  sizeBasisDE ] = zipConformingNew( PHUTelemDE, dimBasisDE, patchBoundaries, p+1, q+1);
    figure
    plotPHTMeshMP(PHUTelem, GIFTmesh)
    
    disp('assign sol Index')
    [PHUTelemDE,solIndexCount] = assignSolIndex2D(PHUTelemDE,patchBoundaries,p+1,q+1);       

    disp('zip conforming')
    [PHUTelem,sizeBasis,numType2Basis,type2Basis] =zipConformingC1_degElev(PHUTelem,PHUTelemDE,GIFTmesh,sizeBasis,patchBoundaries,solIndexCount,p,q);

    %    Material properties
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
    [ stiff, rhs ] = assembleGalerkinSysGIFTMP_degElev( PHUTelem, GIFTmesh, sizeBasis, p, q, Cmat);
    
    disp('Imposing boundary conditions...')
    [ stiff, rhs, bcdof, bcval ] = imposeDirichletNeuPlateHole_degElev(stiff, rhs, PHUTelem, GIFTmesh, p, q, rad, tx);
    
    % pause
    %toc
    disp('Solving the linear system...')
    stepCount
    
    %cond(full(stiff))
    sol0 = stiff\rhs;   
    size(sol0)    
%    pause
    toc
    disp('Estimating the error...')
    [quadRef,estErrorGlobTotal]=recoverDerivEstGalMP2alphaAll3(PHUTelem, GIFTmesh, sol0, target_rel_error, quadList, p, q, Cmat, targetScale);
    estErrorGlobTotal
    
    disp('Computing the errors...')
    [l2relerr, h1relerr]=calcErrorNormsPlateHole_modifiedC_degElev( sol0, PHUTelem, GIFTmesh, p, q, Cmat, Emod, nu, rad, tx );
    %     l2relerr
    %     h1relerr
    
    disp('Plotting the errors...')
    plotErrorPHTElasticPlateHoleC1_degElev( sol0, PHUTelem, GIFTmesh, p, q, Cmat, rad, tx )
    toc
      
    patchA = patchBoundaries{1};
    patchB = patchBoundaries{2};
    edgeA = patchBoundaries{3};
    edgeB = patchBoundaries{4};

    [elementsA] = sortEdgeElem( PHUTelemDE{patchA}, edgeA);
    [elementsB] = sortEdgeElem( PHUTelemDE{patchB}, edgeB);
    
    elementsA=elementsA{1};
    elementsB=elementsB{1};
    
    %%uniform refinement
    for patchIndex=1:numPatches
        quadRef{patchIndex} = 1:size(quadList{patchIndex},1);
        [~, PHUTelemDE{patchIndex}, dimBasisDE(patchIndex)] = refineMesh(quadRef{patchIndex}, quadList{patchIndex}, PHUTelemDE{patchIndex}, p+1, q+1, dimBasisDE(patchIndex));
        [quadList{patchIndex}, PHUTelem{patchIndex}, dimBasis(patchIndex)] = refineMesh(quadRef{patchIndex}, quadList{patchIndex}, PHUTelem{patchIndex}, p, q, dimBasis(patchIndex));    %     end
    end
    
end
