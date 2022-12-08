function [PHUTelemType2BasisInfo] = selectType2BasisAndModifiedC_DE(PHUTelem,type2Basis)

%select type2 basis and modified C from PHUTelem
%use in degree elevation
numPatches=length(PHUTelem);
PHUTelemType2BasisInfo = cell(numPatches, 1);
for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch}) 
        [~,loc]=ismember(type2Basis,PHUTelem{indexPatch}(indexElem).nodesGlobal);  
        numBasis=length(PHUTelem{indexPatch}(indexElem).nodesGlobal);
        tempNodes=setdiff([1:numBasis],loc);
        type2Nodes=PHUTelem{indexPatch}(indexElem).nodesGlobal; 
        type2Nodes(tempNodes)=[];
        PHUTelemType2BasisInfo{indexPatch}(indexElem).type2Nodes=type2Nodes;
        PHUTelemType2BasisInfo{indexPatch}(indexElem).type2NodesC=PHUTelem{indexPatch}(indexElem).modifiedC;
        PHUTelemType2BasisInfo{indexPatch}(indexElem).type2NodesC(tempNodes,:)=[];

    end
end

end

