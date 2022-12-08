function [ PHUTelem ,type2Basis,sizeBasis] = assignNewNodesGlobal( PHUTelem,GIFTmesh,patchBoundaries ,sizeBasis,type2basisNodes,p,q)
right_nodes2 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
right_nodes1 = (p+1):(p+1):(p+1)*(q+1);

left_nodes1 = 1:(p+1):(1+(p+1)*q);
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

down_nodes1=1:(p+1);
down_nodes2=down_nodes1+(p+1);

up_nodes1=(p+1)*(q+1)-p:(p+1)*(q+1);
up_nodes2=up_nodes1-(p+1);

numType2Basis=length(type2basisNodes);

[boundaryInfo,segInfo,subSegInfo,numSeg,patchNeighborInfo,patchEdgeInfo] = createBoundaryInfo3(PHUTelem, patchBoundaries );
removedNodes=[];


for indexPatch=1:length(PHUTelem)
    
    type2Edge=patchEdgeInfo{indexPatch};
    
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

for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        PHUTelem{indexPatch}(indexElem).extraNodes=[];
    end
end

%%%%assign type 2 basis to zeros location of temp nodes global
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
                if length(type2basisNodes)<length(type2Loc)
                    
                    numExtraNodes=length(type2Loc)-length(type2basisNodes);
                    type2Loc(end-numExtraNodes+1:end)=[];
                    tempNodesGlobal=PHUTelem{indexPatch}(indexElem).tempNodesGlobal;
                    tempNodesGlobal(type2Loc)=type2basisNodes;
                    extraNodes=find(tempNodesGlobal==0);
                    tempNodesGlobal(extraNodes)=[];
                    % remove extra loc from tempNodesGlobal
                    PHUTelem{indexPatch}(indexElem).tempNodesGlobal=tempNodesGlobal;
                    PHUTelem{indexPatch}(indexElem).extraNodes=extraNodes;
                else
                    
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
end


removedNodes=unique(removedNodes)
%pause
%re-assign global nodes
temp_nodes=1:sizeBasis+numType2Basis;
replacePattern = temp_nodes;
temp_nodes(removedNodes)=[];
sizeBasis=length(temp_nodes);
replacePattern(temp_nodes) = 1:sizeBasis;

type2Basis=replacePattern(type2basisNodes);
for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(indexElem).children)
            PHUTelem{indexPatch}(indexElem).nodesGlobal = replacePattern(PHUTelem{indexPatch}(indexElem).nodesGlobal);
        end
    end
end

%=================== end assign new nodes global =======================


% disp('after replace patter')
% PHUTelem{1}(1).nodesGlobal
% PHUTelem{1}(2).nodesGlobal
% PHUTelem{1}(3).nodesGlobal
% PHUTelem{1}(4).nodesGlobal
% 
% PHUTelem{2}(1).nodesGlobal
% PHUTelem{2}(2).nodesGlobal
% PHUTelem{2}(3).nodesGlobal
% PHUTelem{2}(4).nodesGlobal
% 
% pause


end

