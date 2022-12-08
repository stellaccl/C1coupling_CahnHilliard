function [PHUTelem,c1SizeBasis ] = removeType2Basis_degElev(PHUTelem,p,q ,sizeBasis)
%remove all type 2 basis functions and c
%for localization method
allNodesGlobal=[];
for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(indexElem).children)
            
            modifiedC=PHUTelem{indexPatch}(indexElem).C;
            tempNodesGlobal=PHUTelem{indexPatch}(indexElem).nodesGlobal;
            PHUTelem{indexPatch}(indexElem).remainNodes=PHUTelem{indexPatch}(indexElem).nodesGlobal;
            
            loc=~ismember(PHUTelem{indexPatch}(indexElem).solIndex,0);
            tempNodesGlobal(loc)=[];
            
            PHUTelem{indexPatch}(indexElem).remainNodes(loc)=0;
            
            modifiedC(loc,:)=[];
            allNodesGlobal=[allNodesGlobal,tempNodesGlobal];
            
            PHUTelem{indexPatch}(indexElem).nodesGlobal=tempNodesGlobal;
            
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
            tempNodesGlobal=PHUTelem{indexPatch}(indexElem).nodesGlobal;
            tempNodesGlobal=replacePattern(tempNodesGlobal);
            PHUTelem{indexPatch}(indexElem).nodesGlobal=tempNodesGlobal;
            
            
            for i=1:length(PHUTelem{indexPatch}(indexElem).remainNodes)
                if PHUTelem{indexPatch}(indexElem).remainNodes(i)~=0
                    
                    PHUTelem{indexPatch}(indexElem).remainNodes(i)=replacePattern(PHUTelem{indexPatch}(indexElem).remainNodes(i));
                    
                end
            end
            
        end
    end
end




end

