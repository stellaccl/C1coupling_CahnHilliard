function [ PHUTelem ,type2basisNodes,sizeBasis] = assignNewNodesGlobal_quadratic( PHUTelem,GIFTmesh,patchBoundaries ,sizeBasis,type2basisNodes,p,q)

[nodes ] = getNodes(p,q);
%[nodesDE ] = getNodes(p+1,q+1);

right_nodes2=nodes.right_nodes2;
right_nodes1=nodes.right_nodes1;

left_nodes1 = nodes.left_nodes1;
left_nodes2 =nodes.left_nodes2;

down_nodes1=nodes.down_nodes1;
down_nodes2=nodes.down_nodes2;

up_nodes1=nodes.up_nodes1;
up_nodes2=nodes.up_nodes2;

%%%%% allocate zero to  type 2 basis location %%%%%%
for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(indexElem).children)
            PHUTelem{indexPatch}(indexElem).tempNodesGlobal=PHUTelem{indexPatch}(indexElem).nodesGlobal;
            loc=~ismember(PHUTelem{indexPatch}(indexElem).solIndex,0);
            PHUTelem{indexPatch}(indexElem).tempNodesGlobal(loc)=0;
            PHUTelem{indexPatch}(indexElem).polyDegree=p*ones(1,size(PHUTelem{indexPatch}(indexElem).nodesGlobal,2));
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
    for indexElem=1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(indexElem).children)
            
            type2Loc=find(PHUTelem{indexPatch}(indexElem).tempNodesGlobal==0);
            if ~isempty(type2Loc)
                if length(type2basisNodes)<length(type2Loc)
                    
                    numExtraNodes=length(type2Loc)-length(type2basisNodes);
                    type2Loc(end-numExtraNodes+1:end)=[];
                    tempNodesGlobal=PHUTelem{indexPatch}(indexElem).tempNodesGlobal;
                    tempNodesGlobal(type2Loc)=type2basisNodes;
                    PHUTelem{indexPatch}(indexElem).polyDegree(type2Loc)=p+1;
                    extraNodes=find(tempNodesGlobal==0);
                    tempNodesGlobal(extraNodes)=[];
                    PHUTelem{indexPatch}(indexElem).polyDegree(extraNodes)=[];
                    % remove extra loc from tempNodesGlobal
                    PHUTelem{indexPatch}(indexElem).tempNodesGlobal=tempNodesGlobal;
                    PHUTelem{indexPatch}(indexElem).extraNodes=extraNodes;
                else
                    
                    additionBasis=length(type2basisNodes)-length(type2Loc);
                    tempNodesGlobal=PHUTelem{indexPatch}(indexElem).tempNodesGlobal;
                    tempNodesGlobal=[tempNodesGlobal,zeros(1,additionBasis)];
                    tempNodesGlobal([type2Loc,up_nodes1(end)+1:up_nodes1(end)+additionBasis])=type2basisNodes;
                    PHUTelem{indexPatch}(indexElem).polyDegree([type2Loc,up_nodes1(end)+1:up_nodes1(end)+additionBasis])=p+1;
                    PHUTelem{indexPatch}(indexElem).tempNodesGlobal=tempNodesGlobal;
                end
            end
            
        end
    end
end






for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(indexElem).children)
            PHUTelem{indexPatch}(indexElem).nodesGlobal=PHUTelem{indexPatch}(indexElem).tempNodesGlobal;
        end
    end
end

% figure
% plotPHTMesh_nodesGlobal2(PHUTelem,GIFTmesh,type2basisNodes,p)
% title('nodesGlobal after assign type 2 basis and zero loc')
% pause

sizeBasis=max(type2basisNodes);

allNodes=[];
for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(indexElem).children)
            allNodes=[allNodes,PHUTelem{indexPatch}(indexElem).nodesGlobal];
        end
    end
end

allNodes=unique(allNodes);
replacePattern=1:sizeBasis;
sizeBasis=length(allNodes);

replacePattern(allNodes)=1:length(allNodes);
type2basisNodes=replacePattern(type2basisNodes);
for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(indexElem).children)
            PHUTelem{indexPatch}(indexElem).nodesGlobal=replacePattern(PHUTelem{indexPatch}(indexElem).nodesGlobal);
            
            for i=1:length(PHUTelem{indexPatch}(indexElem).remainNodes)
                if PHUTelem{indexPatch}(indexElem).remainNodes(i)~=0
                    PHUTelem{indexPatch}(indexElem).remainNodes(i)=replacePattern(PHUTelem{indexPatch}(indexElem).remainNodes(i));
                end
            end
        end
    end
end

% type2basisNodes
%
% figure
% plotPHTMesh_nodesGlobal2(PHUTelem,GIFTmesh,type2basisNodes,p)
% title('nodesGlobal after replace pattern')
% pause


end

