function  [ PHTelem ,type2Basis,sizeBasis] = assignNewNodesGlobal1D( PHTelem,type2Basis,p)
% only for couling between two patch 
%assign zeros loc for type 2 basis
for indexPatch=1:length(PHTelem)
    for indexElem=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(indexElem).children)
            PHTelem{indexPatch}(indexElem).tempNodesGlobal=PHTelem{indexPatch}(indexElem).nodesGlobal;
            loc=find(PHTelem{indexPatch}(indexElem).solIndex~=0);
            if ~isempty(loc)
                PHTelem{indexPatch}(indexElem).tempNodesGlobal(loc)=0;
            end
        end
    end
end

for indexPatch=1:length(PHTelem)
    for indexElem=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(indexElem).children)
            PHTelem{indexPatch}(indexElem).tempNodesGlobal=PHTelem{indexPatch}(indexElem).nodesGlobal;
            PHTelem{indexPatch}(indexElem).polyDegree=p*ones(1,size(PHTelem{indexPatch}(indexElem).nodesGlobal,2));
            loc=find(PHTelem{indexPatch}(indexElem).solIndex~=0);
            if ~isempty(loc)
                PHTelem{indexPatch}(indexElem).tempNodesGlobal(loc)=0;
            end
        end
    end
end

for indexPatch=1:length(PHTelem)
    for indexElem=1:length(PHTelem{indexPatch})
        PHTelem{indexPatch}(indexElem).extraNodes=[];
    end
end

right_nodes=p+1;

for indexPatch=1:length(PHTelem)
    for indexElem=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(indexElem).children)
            
            type2Loc=find(PHTelem{indexPatch}(indexElem).tempNodesGlobal==0);
            if ~isempty(type2Loc)
                if length(type2Basis)<length(type2Loc)
                    
                    numExtraNodes=length(type2Loc)-length(type2Basis);
                    type2Loc(end-numExtraNodes+1:end)=[];
                    tempNodesGlobal=PHTelem{indexPatch}(indexElem).tempNodesGlobal;
                    tempNodesGlobal(type2Loc)=type2Basis;
                    PHTelem{indexPatch}(indexElem).polyDegree(type2Loc)=p+1;
                    extraNodes=find(tempNodesGlobal==0);
                    tempNodesGlobal(extraNodes)=[];
                    PHTelem{indexPatch}(indexElem).polyDegree(extraNodes)=[];
                    % remove extra loc from tempNodesGlobal
                    PHTelem{indexPatch}(indexElem).tempNodesGlobal=tempNodesGlobal;
                    PHTelem{indexPatch}(indexElem).extraNodes=extraNodes;
                else
                    
                    additionBasis=length(type2Basis)-length(type2Loc);
                    tempNodesGlobal=PHTelem{indexPatch}(indexElem).tempNodesGlobal;
                    tempNodesGlobal=[tempNodesGlobal,zeros(1,additionBasis)];
                    tempNodesGlobal([type2Loc,right_nodes(end)+1:right_nodes(end)+additionBasis])=type2Basis;
                    PHTelem{indexPatch}(indexElem).polyDegree([type2Loc,right_nodes(end)+1:right_nodes(end)+additionBasis])=p+1;
                    PHTelem{indexPatch}(indexElem).tempNodesGlobal=tempNodesGlobal;
                end
            end
            
        end
    end
end

for indexPatch=1:length(PHTelem)
    for indexElem=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(indexElem).children)
            PHTelem{indexPatch}(indexElem).nodesGlobal=PHTelem{indexPatch}(indexElem).tempNodesGlobal;
        end
    end
end

sizeBasis=max(type2Basis);

allNodes=[];
for indexPatch=1:length(PHTelem)
    for indexElem=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(indexElem).children)
            allNodes=[allNodes,PHTelem{indexPatch}(indexElem).nodesGlobal];
        end
    end
end

allNodes=unique(allNodes);
replacePattern=1:sizeBasis;
sizeBasis=length(allNodes);

replacePattern(allNodes)=1:length(allNodes);
type2Basis=replacePattern(type2Basis);
for indexPatch=1:length(PHTelem)
    for indexElem=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(indexElem).children)
            PHTelem{indexPatch}(indexElem).nodesGlobal=replacePattern(PHTelem{indexPatch}(indexElem).nodesGlobal);
        end
    end
end



end

