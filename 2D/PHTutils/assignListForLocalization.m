function [ PHUTelem] = assignListForLocalization( PHUTelem,patchBoundaries )

numBoundaries = size(patchBoundaries,1);
for boundaryIndex = 1:numBoundaries
    
    patchAList = patchBoundaries{boundaryIndex,1};
    patchB = patchBoundaries{boundaryIndex,2};
    edgeAList = patchBoundaries{boundaryIndex,3};
    edgeBList = patchBoundaries{boundaryIndex,4};
    
    for indexPatch=1:length(patchAList)
        
        patchA = patchAList(indexPatch);
        edgeA = edgeAList(indexPatch)
        edgeB = edgeBList(indexPatch)
        
      [ elementsA ] = sortEdgeElem(  PHUTelem{patchA}, edgeA)
       [ elementsB ] = sortEdgeElem( PHUTelem{patchB}, edgeB)
        elementsA =cell2mat(elementsA)
        elementsB =cell2mat(elementsB)
%     pause
        
    end
    
end


end
