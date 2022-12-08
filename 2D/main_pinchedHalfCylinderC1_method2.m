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

PHUTelem = cell(numPatches, 1);
for i=1:numPatches
    [PHUTelem{i}, dimBasis(i)] = initPHTmeshGen(p,q,2,2);
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
   
%     figure
%     plotPHTShellMeshMP(PHUTelem, GIFTmesh)
   
    disp('assign sol Index')
    [PHUTelem,solIndexCount] = assignSolIndex2D_quadratic(PHUTelem,patchBoundaries,sizeBasis,p,q);

%     figure
%     plotPHTMesh_solIndexShell( PHUTelem,GIFTmesh,p);
%     title('sol Index')
    %pause

%     figure
%     plotPHTMesh_nodesGlobalShell(PHUTelem,GIFTmesh,p)
%     title('nodes Global')
    %pause
    
    disp('compute matrix m')
    [PHUTelem,m] =zipConformingC1_shell_method2(PHUTelem,GIFTmesh,patchBoundaries,solIndexCount,p,q);


     disp('solving for coef sol')
     [coefSol] = zipConformingC1_boundaryCondition_pinchedCylinder(PHUTelem,m,p,q);
   
   disp('assign new nodes global and modify c')
    [PHUTelem,sizeBasis,type2Basis] = zipConformingC1_quaddratic_modifyC(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,coefSol,solIndexCount,p,q);
  
%     PlotBasisShell( PHUTelem,GIFTmesh,p,q,type2Basis)
    
    
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
    stepCount
    %size(sol)
    vertical_displacement = sol(3*forcedNode-1);
    vertical_displacement_error = vertical_displacement+1.8248e-5;
    norm_vertical_disp_err=norm(vertical_displacement_error)
    
    %sol = zeros(size(sol));
%     figure
      PlotDispShellC1( PHUTelem, GIFTmesh,  p, q, sol)
%     title(['step',num2str(stepCount)])
%     axis equal
%     axis off
%     pause
    % close all
    
    %uniform refinement
    for patchIndex=1:numPatches
        quadRef{patchIndex} = 1:size(quadList{patchIndex},1);
        [quadList{patchIndex}, PHUTelem{patchIndex}, dimBasis(patchIndex)] = refineMesh(quadRef{patchIndex}, quadList{patchIndex}, PHUTelem{patchIndex}, p, q, dimBasis(patchIndex));
    end
end