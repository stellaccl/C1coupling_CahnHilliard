restoredefaultpath
close all
clear all
addpaths

tic;

numPatches = 2;
W = 1;
L = 1;
%GIFTmesh = init2DGeometryGIFTMP('bilinear_rectangle', L, W, numPatches);
%GIFTmesh = init2DGeometryGIFTMP('2patch_leftRight', L, W, numPatches);
GIFTmesh = init2DGeometryGIFTMP('2patch_leftRight', L, W, numPatches);
%pause
p=3;
q=3;

numElemU=2;
numElemV=2;

dimBasis = zeros(1, numPatches);
PHUTelem = cell(numPatches, 1);
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


patchBoundaries = {1,2,2,4};

numSteps =3;
for stepCount = 1:numSteps
    
    [PHUTelem, dimBasis, quadList ] = checkConforming_IGAPack( PHUTelem, dimBasis, patchBoundaries, p, q, quadList );
    [PHUTelem,  sizeBasis ] = zipConformingNew( PHUTelem, dimBasis, patchBoundaries, p, q);
    
    figure
    plotPHTMeshMP(PHUTelem, GIFTmesh)
    %pause
    
    disp('assign sol Index')
    [PHUTelem,solIndexCount] = assignSolIndex2D_quadratic(PHUTelem,patchBoundaries,sizeBasis,p,q);
    
%     figure
%     plotPHTMesh_solIndex(PHUTelem,GIFTmesh,p)
%     title('sol Index')
%     pause
%     %
%     figure
%     plotPHTMesh_nodesGlobal(PHUTelem,GIFTmesh,p)
%     title('nodes Global')
%     pause
    
    disp('compute matrix m')
    [PHUTelem,m] =zipConformingC1_quadratic(PHUTelem,GIFTmesh,patchBoundaries,solIndexCount,p,q);
    %[PHUTelem,m] =zipConformingC1_method3(PHUTelem,GIFTmesh,patchBoundaries,solIndexCount,p,q);
    %
    
    
    disp('apply boundary condition')
    [coefSol] = zipConformingC1_BC_bilinearPatches(PHUTelem,m,p,q );
    
    disp('assign new nodes global and modify c')
    [PHUTelem,sizeBasis,type2Basis] = zipConformingC1_quaddratic_modifyC(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,coefSol,solIndexCount,p,q);
    
    disp('localising type 2 basis function')
    [PHUTelem] = localizeType2Basis2(PHUTelem,type2Basis);
    %testPlotBasisPhys_GIFTmesh3( PHUTelem,GIFTmesh,p,q,3)
    %testPlotBasisPhys_GIFTmesh(PHUTelem,GIFTmesh,p,q,type2Basis)
    
    %testPlotBasisPhys_GIFTmesh(PHUTelem,GIFTmesh,p,q,type2Basis)
    
    %     for i=1:sizeBasis
    %
    %         testPlotBasisPhys_GIFTmesh3(PHUTelem,GIFTmesh,p,q,i)
    %
    %         pause
    %         close all
    %     end
    
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
%    pause
    
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