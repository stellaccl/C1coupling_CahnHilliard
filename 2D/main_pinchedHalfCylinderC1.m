close all
clear all
restoredefaultpath
addpaths

numPatches = 2;
W = 100;
L = 100;

GIFTmesh = init2DGeometryGIFTMP('pinchedHalfCylinder2patch', L, W, numPatches);

dimBasis = zeros(1, numPatches);
p=3;
q=3;
numElemU=2;
numElemV=2;
PHUTelem = cell(numPatches, 1);
for i=1:numPatches
    [PHUTelem{i}, dimBasis(i)] = initPHTmeshGen(p,q,numElemU,numElemV);
    quadList{i} = 1:4;
end

% figure
% plotPHTShellMeshMP(PHUTelem, GIFTmesh)
% axis equal
% pause
patchBoundaries = {1,2,2,4};

numSteps=4;
for stepCount = 1:numSteps
    
    
    [PHUTelem,sizeBasis] = zipConformingNew(PHUTelem,dimBasis,patchBoundaries,p,q);
    [PHUTelem,solIndexCount,patchInfo] = assignSolIndex2D(PHUTelem,patchBoundaries,p,q);
    
%     figure
%     plotPHTShellMeshMP(PHUTelem, GIFTmesh)
    %
    
   
    
%         figure
%         plotPHTMesh_solIndexShell( PHUTelem,GIFTmesh,p);
% %         pause
    
    [boundaryInfo,segInfo,subSegInfo,numSeg,patchNeighborInfo,patchEdgeInfo] = createBoundaryInfo3(PHUTelem, patchBoundaries );
    allSeg=[segInfo;subSegInfo.u;subSegInfo.v];
    
    disp('computing matrix m')
    [m] = computeMatrixM_shell(PHUTelem,GIFTmesh,solIndexCount,allSeg,numSeg,p,q);

   disp('solving for coef sol')
   [coefSol] = zipConformingC1_boundaryCondition_pinchedCylinder(PHUTelem,m,p,q);
    
   disp('modifying c')
   [PHUTelem,sizeBasis,type2Basis] = zipConformingC1_modifyC(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,coefSol,p,q);
    
    %========================= analysis part =============================
    
    % Material properties
    E  = 3e6;
    nu = 0.3;
    t  = 3; % thickness
    
    % constitutive matrix
    memStiff = E*t/(1-nu^2);
    benStiff = E*t^3/12/(1-nu^2);
    
    disp('Assembling the linear system...')
    [ stiff, rhs ] = assembleGalerkin_shellC1( PHUTelem, GIFTmesh, sizeBasis, p, q, nu, memStiff, benStiff );
    condest(stiff)
    
    disp('Imposing boundary conditions...')
    [stiff, rhs, bcdof, bcval, forcedNode ] = imposeDirichlet_pinchedHalfCylinderC1(stiff, rhs, PHUTelem, p, q,sizeBasis,type2Basis);
    condest(stiff)
    
    disp('Solving the linear system...')
    sol = stiff\rhs;

    size(sol)
    vertical_displacement = sol(3*forcedNode-1);
    vertical_displacement_error = vertical_displacement+1.8248e-5;
    norm_vertical_disp_err=norm(vertical_displacement_error)

    figure
    PlotDispShellC1( PHUTelem, GIFTmesh,  p, q, sol)
    title(['step',num2str(stepCount)])
    axis equal
    axis off
%     pause
    % close all
    
    %uniform refinement
    for patchIndex=1:numPatches
        quadRef{patchIndex} = 1:size(quadList{patchIndex},1);
        [quadList{patchIndex}, PHUTelem{patchIndex}, dimBasis(patchIndex)] = refineMesh(quadRef{patchIndex}, quadList{patchIndex}, PHUTelem{patchIndex}, p, q, dimBasis(patchIndex));
    end
    

end