function [PHUTelem,c1SizeBasis ] = removeType2Basis(PHUTelem,p,q ,sizeBasis)
%remove all type 2 basis functions and c
%for localization method
allNodesGlobal=[];
for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(indexElem).children)
            
            modifiedC=PHUTelem{indexPatch}(indexElem).C;
            tempNodesGlobal=PHUTelem{indexPatch}(indexElem).nodesGlobal;
            loc=~ismember(PHUTelem{indexPatch}(indexElem).oriSolIndex,0);
            tempNodesGlobal(loc)=[];
            modifiedC(loc,:)=[];
            allNodesGlobal=[allNodesGlobal,tempNodesGlobal];
            PHUTelem{indexPatch}(indexElem).C1NodesGlobal=tempNodesGlobal;
            PHUTelem{indexPatch}(indexElem).modifiedC=modifiedC;
            PHUTelem{indexPatch}(indexElem).polyDegree=p*ones(1,length(tempNodesGlobal));
            
        end
    end
end
allNodesGlobal=unique(allNodesGlobal);
allNodesGlobal=sort(allNodesGlobal);
% %reassign nodesGlobal
replacePattern=1:sizeBasis;
replacePattern(allNodesGlobal)=1:length(allNodesGlobal);
c1SizeBasis=length(allNodesGlobal);
for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(indexElem).children)
            tempNodesGlobal=PHUTelem{indexPatch}(indexElem).C1NodesGlobal;
            tempNodesGlobal=replacePattern(tempNodesGlobal);
            PHUTelem{indexPatch}(indexElem).C1NodesGlobal=tempNodesGlobal;
        end
    end
end




end

