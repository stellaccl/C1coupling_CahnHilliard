function [ PHUTelem,sizeBasis,numType2Basis,type2Basis] = zipConformingC1_3( PHUTelem,GIFTmesh,solIndexCount,sizeBasis, patchBoundaries,numPatches,numberElementsU,numberElementsV, p, q )
%with compute matrix M 3 for sphere

right_nodes2 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
right_nodes1 = (p+1):(p+1):(p+1)*(q+1);

left_nodes1 = 1:(p+1):(1+(p+1)*q);
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

down_nodes1=1:(p+1);
down_nodes2=down_nodes1+(p+1);

up_nodes1=(p+1)*(q+1)-p:(p+1)*(q+1);
up_nodes2=up_nodes1-(p+1);

[boundaryInfo,segInfo,subSegInfo,numSeg,patchNeighborInfo,patchEdgeInfo] = createBoundaryInfo3(PHUTelem, patchBoundaries );

[m] = computeMatrixM_3(PHUTelem,GIFTmesh,solIndexCount,patchBoundaries,subSegInfo,numSeg,numberElementsU,numberElementsV,p,q );
coefSol = nullMDS(m);

numType2Basis=size(coefSol,2);

type2basisNodes=sizeBasis+1:sizeBasis+numType2Basis;

removedNodes=[];

%=======================assign tempNodesGlobal (zeros at type 2 basis location) ===========================
for indexPatch=1:length(PHUTelem)
    
    type2Edge=patchEdgeInfo{indexPatch};
    
    boundaryElem=boundaryInfo{indexPatch}(:,5);
    boundaryElem=unique(boundaryElem);
    
    for indexElem=1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(indexElem).children)
            
            modifyBasisIndex=[];
            neighbor_left=PHUTelem{indexPatch}(indexElem).neighbor_left;
            neighbor_right=PHUTelem{indexPatch}(indexElem).neighbor_right;
            neighbor_up=PHUTelem{indexPatch}(indexElem).neighbor_up;
            neighbor_down=PHUTelem{indexPatch}(indexElem).neighbor_down;
            for indexEdge=1:length(type2Edge)
                switch type2Edge(indexEdge)
                    case 1
                        if isempty(neighbor_down)
                            modifyBasisIndex=[modifyBasisIndex,down_nodes2,down_nodes1];
                        end
                    case 2
                        if isempty(neighbor_right)
                            modifyBasisIndex=[modifyBasisIndex,right_nodes2,right_nodes1];
                        end
                    case 3
                        if isempty(neighbor_up)
                            modifyBasisIndex=[modifyBasisIndex,up_nodes2,up_nodes1];
                        end
                    case 4
                        if isempty(neighbor_left)
                            modifyBasisIndex=[modifyBasisIndex,left_nodes2,left_nodes1];
                        end
                end
            end
            
            PHUTelem{indexPatch}(indexElem).tempNodesGlobal=PHUTelem{indexPatch}(indexElem).nodesGlobal;
            PHUTelem{indexPatch}(indexElem).tempNodesGlobal( modifyBasisIndex)=0;
            removedNodes=[removedNodes,PHUTelem{indexPatch}(indexElem).nodesGlobal( modifyBasisIndex)];
            
        end
    end
end

disp('finish assign tempNodes Global')

for indexPatch=1:length(PHUTelem)
    
    type2Edge=patchEdgeInfo{indexPatch};
    type2Elem=[];
    for indexEdge=1:length(type2Edge)
        tempElem=sortEdgeElem( PHUTelem{indexPatch}, type2Edge(indexEdge));
        tempElem=cell2mat(tempElem);
        type2Elem = [type2Elem,tempElem];
    end
    
    type2Elem=unique(type2Elem);
    
    for indexElem=1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(indexElem).children)
            if ismember(indexElem,type2Elem)
                
                type2Loc=find(PHUTelem{indexPatch}(indexElem).tempNodesGlobal==0);
                additionBasis=length(type2basisNodes)-length(type2Loc);
                tempNodesGlobal=PHUTelem{indexPatch}(indexElem).tempNodesGlobal;
                tempNodesGlobal=[tempNodesGlobal,zeros(1,additionBasis)];
                tempNodesGlobal([type2Loc,up_nodes1(end)+1:up_nodes1(end)+additionBasis])=type2basisNodes;
                PHUTelem{indexPatch}(indexElem).tempNodesGlobal=tempNodesGlobal;
            end
            
            PHUTelem{indexPatch}(indexElem).nodesGlobal=PHUTelem{indexPatch}(indexElem).tempNodesGlobal;
        end
    end
    
end

removedNodes=unique(removedNodes);
%re-assign global nodes
temp_nodes=1:sizeBasis+numType2Basis;
replacePattern = temp_nodes;
temp_nodes(removedNodes)=[];
sizeBasis=length(temp_nodes);
replacePattern(temp_nodes) = 1:sizeBasis;

type2Basis=replacePattern(type2basisNodes);

for indexPatch=1:numPatches
    for indexElem=1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(indexElem).children)
            PHUTelem{indexPatch}(indexElem).nodesGlobal = replacePattern(PHUTelem{indexPatch}(indexElem).nodesGlobal);
        end
    end
end

%=================== end assign new nodes global =======================

[PHUTelem] = modifyC_3( PHUTelem,coefSol,numType2Basis,patchEdgeInfo,boundaryInfo,type2Basis,p,q);
end

