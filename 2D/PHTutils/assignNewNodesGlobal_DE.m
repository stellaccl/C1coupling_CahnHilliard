function [PHUTelem] = assignNewNodesGlobal_DE(PHUTelem,PHUTelemType2BasisInfo,PHUTelemAdditionalBasisInfo,GIFTmesh,patchBoundaries ,sizeBasis,type2basisNodes,p,q)


%change size of modified C such that match with the Degree elevated size
for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        tempC=PHUTelem{indexPatch}(indexElem).modifiedC;
        newC=zeros(size(tempC,1),(p+2)*(q+2));
        PHUTelem{indexPatch}(indexElem).modifiedC=newC;
        PHUTelem{indexPatch}(indexElem).modifiedC(:,1:size(tempC,2))=tempC;
    end
end


%add type 2 basis and type2NodesC
%add
for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        
        type2Nodes=PHUTelemType2BasisInfo{indexPatch}(indexElem).type2Nodes;
        type2NodesC=PHUTelemType2BasisInfo{indexPatch}(indexElem).type2NodesC;
        
        additionalBasis=PHUTelemAdditionalBasisInfo{indexPatch}(indexElem).additionalBasis;
        additionalBasisC=PHUTelemAdditionalBasisInfo{indexPatch}(indexElem).additionalBasisC;
        
        numBasis=length(PHUTelem{indexPatch}(indexElem).nodesGlobal);
        PHUTelem{indexPatch}(indexElem).polyDegree=p*ones(1,numBasis);
        
        if ~isempty(type2Nodes)
            numType2Basis=length(type2Nodes);
            numBasis=length(PHUTelem{indexPatch}(indexElem).nodesGlobal);
            PHUTelem{indexPatch}(indexElem).nodesGlobal(numBasis+1:numBasis+numType2Basis)=type2Nodes;
            PHUTelem{indexPatch}(indexElem).modifiedC(numBasis+1:numBasis+numType2Basis,:)=type2NodesC;
            PHUTelem{indexPatch}(indexElem).polyDegree(numBasis+1:numBasis+numType2Basis)=p+1;
        end
        
        if ~isempty(additionalBasis)
            numAdditionalBasis=length(additionalBasis);
            numBasis=length(PHUTelem{indexPatch}(indexElem).nodesGlobal);
            PHUTelem{indexPatch}(indexElem).nodesGlobal(numBasis+1:numBasis+numAdditionalBasis)=additionalBasis;
            PHUTelem{indexPatch}(indexElem).modifiedC(numBasis+1:numBasis+numAdditionalBasis,:)=additionalBasisC;
            PHUTelem{indexPatch}(indexElem).polyDegree(numBasis+1:numBasis+numAdditionalBasis)=p+1;
        end

    end
end

end

