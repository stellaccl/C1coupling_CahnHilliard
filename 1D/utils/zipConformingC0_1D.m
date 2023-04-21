function [ PHUTelem, sizeBasis ] =zipConformingC0_1D( PHUTelem )
%only for 2 patch

indexPatch=1;
allNodes=[];
for indexElem=1:length(PHUTelem{indexPatch})
    if isempty(PHUTelem{indexPatch}(indexElem).children)
        PHUTelem{indexPatch}(indexElem).nodesGlobal=PHUTelem{indexPatch}(indexElem).nodes;
        allNodes=[allNodes,PHUTelem{indexPatch}(indexElem).nodesGlobal];
    end
end

count=max(allNodes);
adj=count-1;
indexPatch=2;
for indexElem=1:length(PHUTelem{indexPatch})
    if isempty(PHUTelem{indexPatch}(indexElem).children)
        PHUTelem{indexPatch}(indexElem).nodesGlobal=PHUTelem{indexPatch}(indexElem).nodes+adj;
         allNodes=[allNodes,PHUTelem{indexPatch}(indexElem).nodesGlobal];
    end
end
sizeBasis =max(allNodes);
end

