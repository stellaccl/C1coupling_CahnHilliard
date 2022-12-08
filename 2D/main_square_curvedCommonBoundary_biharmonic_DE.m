restoredefaultpath
close all
clear all
addpaths

tic;

numPatches = 2;
W = 1;
L = 1;
%GIFTmesh = init2DGeometryGIFTMP('rectangle', L, W, numPatches);
GIFTmesh = init2DGeometryGIFTMP('2patch_curvedCommonBoundary', L, W, numPatches);


%%%%%%%%%%%%%%%%% initialize %%%%%%%%%%%%%%%%%%
p=3;
q=3;

numElemU=3;
numElemV=3;

dimBasis = zeros(1, numPatches);
PHUTelem = cell(numPatches, 1);
quadList = cell(numPatches, 1);

dimBasisDE = zeros(1, numPatches);
PHUTelemDE = cell(numPatches, 1);
quadListDE = cell(numPatches, 1);

if p==2
    
    for indexPatch=1:numPatches
        [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmesh_quadratic( p,q,numElemU,numElemV);
        numElem=length(PHUTelem{indexPatch});
        quadList{indexPatch} =1:numElem;
        
    end
    
    for indexPatch=1:numPatches
        [PHUTelemDE{indexPatch}, dimBasisDE(indexPatch)] = initPHTmeshGen( p+1,q+1,numElemU,numElemV);
        numElem=length(PHUTelemDE{indexPatch});
        quadListDE{indexPatch} =1:numElem;
    end
    
else
    
    for indexPatch=1:numPatches
        [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmeshGen( p,q,numElemU,numElemV);
        numElem=length(PHUTelem{indexPatch});
        quadList{indexPatch} =1:numElem;
        
        [PHUTelemDE{indexPatch}, dimBasisDE(indexPatch)] = initPHTmeshGen( p+1,q+1,numElemU,numElemV);
        numElem=length(PHUTelemDE{indexPatch});
        quadListDE{indexPatch} =1:numElem;
    end
    
end



%     figure
%     plotPHTMeshMP(PHUTelem, GIFTmesh)
%    pause

%%%%%%%%%%%%%%%%%%%%%%% zipConforming

patchBoundaries = {1,2,2,4};

numSteps =4;
for stepCount = 1:numSteps
    
    
   % figure
   % plotPHTMeshMP(PHUTelem, GIFTmesh)
    %  pause
    
    [PHUTelem, dimBasis, quadList ] = checkConforming_IGAPack( PHUTelem, dimBasis, patchBoundaries, p, q, quadList );
    [PHUTelem,  sizeBasis ] = zipConformingNew( PHUTelem, dimBasis, patchBoundaries, p, q);
    
    [PHUTelemDE, dimBasisDE, quadListDE ] = checkConforming_IGAPack( PHUTelemDE, dimBasisDE, patchBoundaries, p+1, q+1, quadListDE );
    [PHUTelemDE,  sizeBasisDE ] = zipConformingNew( PHUTelemDE, dimBasisDE, patchBoundaries, p+1, q+1);
    
    disp('assign sol Index')
    [PHUTelem,~] = assignSolIndex2D_DE_2(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,p,q);
    [PHUTelemDE,solIndexCount] = assignSolIndex2D_DE_2(PHUTelemDE,GIFTmesh,patchBoundaries,sizeBasisDE,p+1,q+1);
    
%          figure
%         plotPHTMesh_solIndex( PHUTelem,GIFTmesh,p)
%         title('sol Index')
%         pause
%     
%     
%          figure
%         plotPHTMesh_solIndex( PHUTelemDE,GIFTmesh,p+1)
%         title('sol Index DE')
%         pause


        
    
    disp('remove type2 basis function')
    [PHUTelem,sizeBasis] = removeType2Basis_degElev(PHUTelem,p,q ,sizeBasis);
    
    
    disp('compute matrix m')
    [PHUTelemDE,m] =zipConformingC1_quadratic(PHUTelemDE,GIFTmesh,patchBoundaries,solIndexCount,p+1,q+1);
    
    disp('apply boundary condition')
    [coefSol] = zipConformingC1_p2_BC_2patchSquareBiharmonic(PHUTelemDE,m,p+1,q+1);

    disp('assign type 2 basis')
    [PHUTelem,sizeBasis,type2Basis] = assignC1basis_degElev(PHUTelem,PHUTelemDE,sizeBasis,coefSol,solIndexCount,p,q);
     
    disp('assign additional basis')
    [PHUTelem,PHUTelemDE,sizeBasis ] = assignAdditionalNodes_degElev( PHUTelem,PHUTelemDE,sizeBasis,p,q);
   
   
    
    Emod=1e5;
    nu=0;
    Cmat = Emod/((1+nu)*(1-2*nu))*[1-nu, nu, 0; nu, 1-nu, 0; 0, 0, (1-2*nu)/2];
    bound_disp = 0.1;
    
    disp('Assembling the linear system...')
    [stiff,rhs] = assembleBiharmonic_c1_degElev(PHUTelem,GIFTmesh,sizeBasis,p,q);
    disp('Imposing boundary conditions...')
    [ stiff, rhs, bcdof, bcval] = imposeDirichletBiharmonic_c1_degElev(stiff, rhs, PHUTelem,PHUTelemDE, p, q ,type2Basis);
    
    
    disp('Solving the linear system...')
    sol0 = stiff\rhs;
    
    figure
    PlotDispPlate_c1_degElev(PHUTelem, GIFTmesh,  p, q, sol0);
    axis equal
 
   size(sol0)
   disp('cal error...')
   [l2relerr] =calcErrorNormsBiharmonic_square_degElev(PHUTelem, GIFTmesh,  p, q, sol0)
   
   
   %%%%%%%%%%%%%%%%%%%%% refinement %%%%%%%%%%%%%%%%%%%%
    if stepCount<numSteps
        numElemU=numElemU*2;
        numElemV=numElemV*2;
        for indexPatch=1:numPatches
            if p==2
                [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmesh_quadratic( p,q,numElemU,numElemV);
                [PHUTelemDE{indexPatch}, dimBasisDE(indexPatch)] = initPHTmeshGen( p+1,q+1,numElemU,numElemV);
            else
                [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmeshGen( p,q,numElemU,numElemV);
                [PHUTelemDE{indexPatch}, dimBasisDE(indexPatch)] = initPHTmeshGen( p+1,q+1,numElemU,numElemV);
            end
            numElem=length(PHUTelem{1});
            quadList{indexPatch} =1:numElem;
            numElem=length(PHUTelemDE{1});
            quadListDE{indexPatch} =1:numElem;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

