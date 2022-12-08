close all
clear all
restoredefaultpath

addpaths
numPatches = 3;
numElemU=4;
numElemV=4;
p=3;
q=3;

GIFTmesh = init2DGeometryGIFTMP('threePatch');

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

patchBoundaries = {1,2,3,1;...
    [1,2],3,[4,4],[1,2]};

numSteps =2;
for stepCount = 1:numSteps

    [PHUTelem, dimBasis, quadList ] = checkConforming_IGAPack( PHUTelem, dimBasis, patchBoundaries, p, q, quadList );
    [PHUTelem,  sizeBasis ] = zipConformingNew( PHUTelem, dimBasis, patchBoundaries, p, q);
    
%     figure
%     plotPHTMeshMP(PHUTelem, GIFTmesh)

    
    % figure
    % plotPHTMesh_nodesGlobal( PHUTelem,GIFTmesh,p)
    % title('nodes Global')
    % pause
    
    disp('assign sol Index')
    [PHUTelem,solIndexCount] = assignSolIndex2D_quadratic(PHUTelem,patchBoundaries,sizeBasis,p,q);
    
    % figure
    % plotPHTMesh_solIndex(PHUTelem,GIFTmesh,p)
    % title('sol Index')
    % pause
    
    % figure
    % plotPHTMesh_solIndex( PHUTelem,GIFTmesh,p)
    % title('sol Index after assign neighbor')
    % pause
    
    disp('compute matrix m')
    [PHUTelem,m] =zipConformingC1_quadratic(PHUTelem,GIFTmesh,patchBoundaries,solIndexCount,p,q);
    disp('apply boundary condition')
    [coefSol] = zipConformingC1_quadratic_boundaryCondition_threePatch(PHUTelem,m,p,q);
    disp('assign new nodes global and modify c')
    [PHUTelem,sizeBasis,type2Basis] = zipConformingC1_quaddratic_modifyC(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,coefSol,solIndexCount,p,q);
    
    %[PHUTelem,sizeBasis,numType2Basis,type2Basis] =zipConformingC1_squareBiharmonic_quartic(PHUTelem,GIFTmesh,sizeBasis,patchBoundaries,solIndexCount,p,q);
    
    disp('localising type 2 basis function')
    [PHUTelem] = localizeType2Basis2(PHUTelem,type2Basis);
    
    % for i=1:sizeBasis
    %
    %     testPlotBasisPhys_GIFTmesh3(PHUTelem,GIFTmesh,p,q,i)
    %
    %     pause
    %     close all
    % end
    
    %disp('localising type 2 basis function')
    %[PHUTelem] = localizeType2Basis2(PHUTelem,type2Basis);
    
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