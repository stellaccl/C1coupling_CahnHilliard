function [m ] = computeMatrixM_3(PHUTelem,GIFTmesh,solIndexCount,patchBoundaries,subSegInfo,numSeg,numberElementsU,numberElementsV,p,q )
%for five patch example
[PHUTelemV,GIFTmeshV,PHUTelemCR,GIFTmeshCR,PHUTelemCL,GIFTmeshCL] = rotateGiftMeshAndPHUTelem2(numberElementsU,numberElementsV);
[PHUTelemTBleft,GIFTmeshTBleft,PHUTelemTBright,GIFTmeshTBright] = rotateGiftMeshAndPHUTelemTopBottom(numberElementsU,numberElementsV);
%[boundaryInfo,segInfo,subSegInfo,numSeg,patchNeighborInfo,patchEdgeInfo] = createBoundaryInfo3(PHUTelem, patchBoundaries );

right_nodes2 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
right_nodes1 = (p+1):(p+1):(p+1)*(q+1);

left_nodes1 = 1:(p+1):(1+(p+1)*q);
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

down_nodes1=1:(p+1);
down_nodes2=down_nodes1+(p+1);

up_nodes1=(p+1)*(q+1)-p:(p+1)*(q+1);
up_nodes2=up_nodes1-(p+1);

numColumns=solIndexCount;
numPts=ceil(numColumns/numSeg);
if numPts<4
    numPts=4;
end
fudge=0;
vref = linspace(-1+fudge,1-fudge,numPts);
uref = linspace(-1+fudge,1-fudge,numPts);

m=zeros(numSeg*numPts,numColumns);

i=0;
for indexBoundary=1:size(patchBoundaries,1)
    
    patch1List = patchBoundaries{indexBoundary,1};
    patch2 = patchBoundaries{indexBoundary,2};
    edge1List = patchBoundaries{indexBoundary,3};
    edge2List = patchBoundaries{indexBoundary,4};
    
    for indexPatch=1:length(patch1List)
        
        patch1 = patch1List(indexPatch);
        edge1 = edge1List(indexPatch);
        edge2 = edge2List(indexPatch);
        
        tempElem1 = sortEdgeElem( PHUTelem{patch1}, edge1);
        tempElem2 = sortEdgeElem( PHUTelem{patch2}, edge2);
        
        tempElem1=cell2mat(tempElem1);
        tempElem2=cell2mat(tempElem2);
        
        for indexSeg=1:length(tempElem1)
            
            elem1=tempElem1(indexSeg);
            elem2=tempElem2(indexSeg);
            
            switch edge1
                case 1
                    solIndexMidP1=PHUTelem{patch1}(elem1).solIndex(down_nodes1);
                    solIndexP1=PHUTelem{patch1}(elem1).solIndex(down_nodes2);
                case 2
                    solIndexMidP1=PHUTelem{patch1}(elem1).solIndex(right_nodes1);
                    solIndexP1=PHUTelem{patch1}(elem1).solIndex(right_nodes2);
                case 3
                    solIndexMidP1=PHUTelem{patch1}(elem1).solIndex(up_nodes1);
                    solIndexP1=PHUTelem{patch1}(elem1).solIndex(up_nodes2);
                case 4
                    solIndexMidP1=PHUTelem{patch1}(elem1).solIndex(left_nodes1);
                    solIndexP1=PHUTelem{patch1}(elem1).solIndex(left_nodes2);
            end
            
            switch edge2
                case 1
                    solIndexMidP2=PHUTelem{patch2}(elem2).solIndex(down_nodes1);
                    solIndexP2=PHUTelem{patch2}(elem2).solIndex(down_nodes2);
                case 2
                    solIndexMidP2=PHUTelem{patch2}(elem2).solIndex(right_nodes1);
                    solIndexP2=PHUTelem{patch2}(elem2).solIndex(right_nodes2);
                case 3
                    solIndexMidP2=PHUTelem{patch2}(elem2).solIndex(up_nodes1);
                    solIndexP2=PHUTelem{patch2}(elem2).solIndex(up_nodes2);
                case 4
                    solIndexMidP2=PHUTelem{patch2}(elem2).solIndex(left_nodes1);
                    solIndexP2=PHUTelem{patch2}(elem2).solIndex(left_nodes2);
            end
            
            if solIndexMidP1 ~= solIndexMidP2
                disp('error')
                pause
            end
            
            if all(ismember([patch1,patch2],[1,2])) || all(ismember([patch1,patch2],[2,3])) || all(ismember([patch1,patch2],[1,4])) || all(ismember([patch1,patch2],[3,4]))
                %v shape
                tempPHUTelemRotated=PHUTelemV;
                tempGIFTmeshRotated=GIFTmeshV;
                
%                 figure
%                 plotPHTMeshMP(tempPHUTelemRotated, tempGIFTmeshRotated);
%                 pause
%                 
                switchSolIndex=0;
            elseif all(ismember([patch1,patch2],[1,5])) || all(ismember([patch1,patch2],[4,5]))
                % center left shape
                tempPHUTelemRotated=PHUTelemCL;
                tempGIFTmeshRotated=GIFTmeshCL;
                switchSolIndex=0;
%                 disp('CL')
%                 patch1
%                 patch2
%                 
%                 figure
%                 plotPHTMeshMP(tempPHUTelemRotated, tempGIFTmeshRotated);
%                 pause
                
            elseif all(ismember([patch1,patch2],[2,5])) || all(ismember([patch1,patch2],[3,5]))
                %center right shape
                tempPHUTelemRotated=PHUTelemCR;
                tempGIFTmeshRotated=GIFTmeshCR;
                switchSolIndex=1;
                
                
%                 disp('CR')
%                 patch1
%                 patch2
%                 
%                 
%                 figure
%                 plotPHTMeshMP(tempPHUTelemRotated, tempGIFTmeshRotated);
%                 pause
                
            elseif all(ismember([patch1,patch2],[1,6])) || all(ismember([patch1,patch2],[4,6]))
                %center right shape
                tempPHUTelemRotated=PHUTelemTBleft;
                tempGIFTmeshRotated=GIFTmeshTBleft;
                switchSolIndex=0;
                
                
%                 disp('TBL')
%                 patch1
%                 patch2
%                 
%                 
%                 figure
%                 plotPHTMeshMP(tempPHUTelemRotated, tempGIFTmeshRotated);
%                 pause
                
            elseif all(ismember([patch1,patch2],[2,6])) || all(ismember([patch1,patch2],[3,6]))
                %center right shape
                tempPHUTelemRotated=PHUTelemTBright;
                tempGIFTmeshRotated=GIFTmeshTBright;
                switchSolIndex=1;
                
%                 disp('TBR')
%                 patch1
%                 patch2
%                 
%                 
%                 figure
%                 plotPHTMeshMP(tempPHUTelemRotated, tempGIFTmeshRotated);
%                 pause
                
            else
                disp('error unable to determine rotated PHUTelem and GIFTmesh')
                pause
            end
            
            elemRotatedPatch1 = sortEdgeElem( tempPHUTelemRotated{1}, 2);
            elemRotatedPatch2 = sortEdgeElem( tempPHUTelemRotated{2}, 4);
            
            elemRotatedPatch1=cell2mat(elemRotatedPatch1);
            elemRotatedPatch2=cell2mat(elemRotatedPatch2);
            
            PHUTelemRotated1=tempPHUTelemRotated{1}(elemRotatedPatch1(indexSeg));
            PHUTelemRotated2=tempPHUTelemRotated{2}(elemRotatedPatch2(indexSeg));
            
            GIFTmeshRotated1=tempGIFTmeshRotated{1};
            GIFTmeshRotated2=tempGIFTmeshRotated{2};
            
            urefP1=1;
            urefP2=-1;
            basisIndexP1=[3 4;7 8;11 12; 15 16];
            basisIndexP2=[1 2;5 6;9 10; 13 14];
            
            if switchSolIndex
                tempSolIndex=[solIndexP2;solIndexMidP1;solIndexP1];
            else
                tempSolIndex=[solIndexP1;solIndexMidP1;solIndexP2];
            end
            
            % tempSolIndex
            % pause
            
            for ic=1:numPts
                i = i+1;
                
                [a,b,c] = computeFunctionABC_GIFTmesh_vDirection2(PHUTelemRotated1,PHUTelemRotated2, GIFTmeshRotated1, GIFTmeshRotated2,vref(ic),urefP1,urefP2);
                a=-a;
                
                [ ~, dRdu_patch2,  dRdv_patch2] = tensorProductBernsteinPolynomial(urefP2,vref(ic),p,q );
                [ ~, dRdu_patch1,  dRdv_patch1] = tensorProductBernsteinPolynomial(urefP1,vref(ic),p,q );
                
                solIndex= tempSolIndex(:,1);
                m(i,solIndex(1))=b*dRdv_patch1(basisIndexP1(1,1))+c*dRdu_patch1(basisIndexP1(1,1));
                m(i,solIndex(2))=b*dRdv_patch1(basisIndexP1(1,2))+c*dRdu_patch1(basisIndexP1(1,2))+a*dRdu_patch2(basisIndexP2(1,1));
                m(i,solIndex(3))=a*dRdu_patch2(basisIndexP2(1,2));
                
                solIndex= tempSolIndex(:,2);
                m(i,solIndex(1))=b*dRdv_patch1(basisIndexP1(2,1))+c*dRdu_patch1(basisIndexP1(2,1));
                m(i,solIndex(2))=b*dRdv_patch1(basisIndexP1(2,2))+c*dRdu_patch1(basisIndexP1(2,2))+a*dRdu_patch2(basisIndexP2(2,1));
                m(i,solIndex(3))=a*dRdu_patch2(basisIndexP2(2,2));
                
                solIndex= tempSolIndex(:,3);
                m(i,solIndex(1))=b*dRdv_patch1(basisIndexP1(3,1))+c*dRdu_patch1(basisIndexP1(3,1));
                m(i,solIndex(2))=b*dRdv_patch1(basisIndexP1(3,2))+c*dRdu_patch1(basisIndexP1(3,2))+a*dRdu_patch2(basisIndexP2(3,1));
                m(i,solIndex(3))=a*dRdu_patch2(basisIndexP2(3,2));
                
                solIndex= tempSolIndex(:,4);
                m(i,solIndex(1))=b*dRdv_patch1(basisIndexP1(4,1))+c*dRdu_patch1(basisIndexP1(4,1));
                m(i,solIndex(2))=b*dRdv_patch1(basisIndexP1(4,2))+c*dRdu_patch1(basisIndexP1(4,2))+a*dRdu_patch2(basisIndexP2(4,1));
                m(i,solIndex(3))=a*dRdu_patch2(basisIndexP2(4,2));
                
            end
            
        end
    end
end
%%=========================================================================
%%=========================== deal with subBrance =========================

subSegV=subSegInfo.v; % sub brance v in u direction
subSegU=subSegInfo.u; %sub brance u in v direction

for indexSeg=1:size(subSegV,1)
    
    patch1=subSegV(indexSeg,1);
    patch2=subSegV(indexSeg,2);
    edge1=subSegV(indexSeg,3);
    edge2=subSegV(indexSeg,4);
    elem1=subSegV(indexSeg,5);
    elem2=subSegV(indexSeg,6);
    
    
    switch edge1
        case 1
            vrefPatch1=-1;
        case 3
            vrefPatch1=1;
            vrefPatch2=-1;
    end
    
    
    if vrefPatch1==-1
        %disp('swith patch')
        patch1=subSegV(indexSeg,2);
        patch2=subSegV(indexSeg,1);
        edge1=subSegV(indexSeg,4);
        edge2=subSegV(indexSeg,3);
        elem1=subSegV(indexSeg,6);
        elem2=subSegV(indexSeg,5);
        
        vrefPatch1=1;
        vrefPatch2=-1;
    end
    
    switch edge1
        case 1
            solIndexMidP1=PHUTelem{patch1}(elem1).solIndex(down_nodes1);
            solIndexP1=PHUTelem{patch1}(elem1).solIndex(down_nodes2);
            basisIndexP1=[down_nodes2;down_nodes1];
        case 2
            solIndexMidP1=PHUTelem{patch1}(elem1).solIndex(right_nodes1);
            solIndexP1=PHUTelem{patch1}(elem1).solIndex(right_nodes2);
            basisIndexP1=[right_nodes2;right_nodes1];
        case 3
            solIndexMidP1=PHUTelem{patch1}(elem1).solIndex(up_nodes1);
            solIndexP1=PHUTelem{patch1}(elem1).solIndex(up_nodes2);
            basisIndexP1=[up_nodes2;up_nodes1];
        case 4
            solIndexMidP1=PHUTelem{patch1}(elem1).solIndex(left_nodes1);
            solIndexP1=PHUTelem{patch1}(elem1).solIndex(left_nodes2);
            basisIndexP1=[left_nodes2;left_nodes1];
    end
    
    switch edge2
        case 1
            solIndexMidP2=PHUTelem{patch2}(elem2).solIndex(down_nodes1);
            solIndexP2=PHUTelem{patch2}(elem2).solIndex(down_nodes2);
            basisIndexP2=[down_nodes1;down_nodes2];
        case 2
            solIndexMidP2=PHUTelem{patch2}(elem2).solIndex(right_nodes1);
            solIndexP2=PHUTelem{patch2}(elem2).solIndex(right_nodes2);
            basisIndexP2=[right_nodes1;right_nodes2];
        case 3
            solIndexMidP2=PHUTelem{patch2}(elem2).solIndex(up_nodes1);
            solIndexP2=PHUTelem{patch2}(elem2).solIndex(up_nodes2);
            basisIndexP2=[up_nodes1;up_nodes2];
        case 4
            solIndexMidP2=PHUTelem{patch2}(elem2).solIndex(left_nodes1);
            solIndexP2=PHUTelem{patch2}(elem2).solIndex(left_nodes2);
            basisIndexP2=[left_nodes1;left_nodes2];
    end
    
    if solIndexMidP1 ~= solIndexMidP2
        disp('error')
        pause
    end
    
    solIndex=[solIndexP1;solIndexMidP1;solIndexP2];
    basisIndex=[basisIndexP1;basisIndexP2];
    for ic=1:numPts
        
        [ ~, dRdu_patch1,  dRdv_patch1] = tensorProductBernsteinPolynomial(uref(ic),vrefPatch1,p,q );
        [ ~, dRdu_patch2,  dRdv_patch2] = tensorProductBernsteinPolynomial(uref(ic),vrefPatch2,p,q );
        
        [dRdu_patch1,dRdv_patch1,dRdu_patch2,dRdv_patch2] =scalingDerivatives2D(PHUTelem,patch1,patch2,elem1,elem2,dRdu_patch1,dRdv_patch1,dRdu_patch2,dRdv_patch2);
        
        [a1,b1,c1] = computeFunctionABC_GIFTmesh_uDirection( PHUTelem,GIFTmesh,patch1,elem1,patch2,elem2,uref(ic));
        a1=-a1;
        
        i = i+1;
        [ m ] = ComputeMatrixM_uDirection2( m, solIndex,basisIndex,i,a1,b1,c1,dRdu_patch1,dRdv_patch1,dRdu_patch2,dRdv_patch2,p,q );
        
    end
end


for indexSeg=1:size(subSegU,1)
    
    patch1=subSegU(indexSeg,1);
    patch2=subSegU(indexSeg,2);
    edge1=subSegU(indexSeg,3);
    edge2=subSegU(indexSeg,4);
    elem1=subSegU(indexSeg,5);
    elem2=subSegU(indexSeg,6);
    
    switch edge1
        case 2
            urefPatch1=1;
        case 4
            urefPatch1=-1;
            urefPatch2=1;
    end
    
    
    if urefPatch1==1
        %disp('swith patch')
        patch1=subSegU(indexSeg,2);
        patch2=subSegU(indexSeg,1);
        edge1=subSegU(indexSeg,4);
        edge2=subSegU(indexSeg,3);
        elem1=subSegU(indexSeg,6);
        elem2=subSegU(indexSeg,5);
        
        urefPatch1=-1;
        urefPatch2=1;
        
    end
    
    
    switch edge1
        case 1
            solIndexMidP1=PHUTelem{patch1}(elem1).solIndex(down_nodes1);
            solIndexP1=PHUTelem{patch1}(elem1).solIndex(down_nodes2);
            basisIndexP1=[down_nodes2;down_nodes1];
        case 2
            solIndexMidP1=PHUTelem{patch1}(elem1).solIndex(right_nodes1);
            solIndexP1=PHUTelem{patch1}(elem1).solIndex(right_nodes2);
            basisIndexP1=[right_nodes2;right_nodes1];
        case 3
            solIndexMidP1=PHUTelem{patch1}(elem1).solIndex(up_nodes1);
            solIndexP1=PHUTelem{patch1}(elem1).solIndex(up_nodes2);
            basisIndexP1=[up_nodes2;up_nodes1];
        case 4
            solIndexMidP1=PHUTelem{patch1}(elem1).solIndex(left_nodes1);
            solIndexP1=PHUTelem{patch1}(elem1).solIndex(left_nodes2);
            basisIndexP1=[left_nodes2;left_nodes1];
    end
    
    switch edge2
        case 1
            solIndexMidP2=PHUTelem{patch2}(elem2).solIndex(down_nodes1);
            solIndexP2=PHUTelem{patch2}(elem2).solIndex(down_nodes2);
            basisIndexP2=[down_nodes1;down_nodes2];
        case 2
            solIndexMidP2=PHUTelem{patch2}(elem2).solIndex(right_nodes1);
            solIndexP2=PHUTelem{patch2}(elem2).solIndex(right_nodes2);
            basisIndexP2=[right_nodes1;right_nodes2];
        case 3
            solIndexMidP2=PHUTelem{patch2}(elem2).solIndex(up_nodes1);
            solIndexP2=PHUTelem{patch2}(elem2).solIndex(up_nodes2);
            basisIndexP2=[up_nodes1;up_nodes2];
        case 4
            solIndexMidP2=PHUTelem{patch2}(elem2).solIndex(left_nodes1);
            solIndexP2=PHUTelem{patch2}(elem2).solIndex(left_nodes2);
            basisIndexP2=[left_nodes1;left_nodes2];
    end
    
    if solIndexMidP1 ~= solIndexMidP2
        disp('error')
        pause
    end
    
    solIndex=[solIndexP1;solIndexMidP1;solIndexP2];
    basisIndex=[basisIndexP1;basisIndexP2];
    for ic=1:numPts
        [ ~, dRdu_patch1,  dRdv_patch1] = tensorProductBernsteinPolynomial(urefPatch1,vref(ic),p,q );
        [ ~, dRdu_patch2,  dRdv_patch2] = tensorProductBernsteinPolynomial(urefPatch2,vref(ic),p,q );
        
        [dRdu_patch1,dRdv_patch1,dRdu_patch2,dRdv_patch2] =scalingDerivatives2D(PHUTelem,patch1,patch2,elem1,elem2,dRdu_patch1,dRdv_patch1,dRdu_patch2,dRdv_patch2);
        [a1,b1,c1] = computeFunctionABC_GIFTmesh_vDirection(  PHUTelem,GIFTmesh,patch1,elem1,patch2,elem2,vref(ic));
        a1=-a1;
        
        i = i+1;
        
        [ m ] = ComputeMatrixM_vDirection2( m, solIndex,basisIndex,i,a1,b1,c1,dRdu_patch1,dRdv_patch1,dRdu_patch2,dRdv_patch2,p,q );
        
    end
end

end

