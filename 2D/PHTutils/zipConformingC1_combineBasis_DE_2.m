function[PHUTelem,PHUTelemDE,sizeBasis,type2Basis,additionalBasis] = zipConformingC1_combineBasis_DE_2(PHUTelem,PHUTelemDE,PHUTelemType2BasisInfo,GIFTmesh,patchBoundaries,sizeBasis,tempType2Basis,p,q)

%initialize assign (modifiedC and polyDegree)
for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        PHUTelem{indexPatch}(indexElem).modifiedC=PHUTelem{indexPatch}(indexElem).C;
        numBasis=length(PHUTelem{indexPatch}(indexElem).nodesGlobal);
        PHUTelem{indexPatch}(indexElem).polyDegree=p*ones(1,numBasis);
    end
end

allBasis=[];
%remove basis from PHUTelem
for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        
        removeNodes=PHUTelem{indexPatch}(indexElem).removeNodes;
        PHUTelem{indexPatch}(indexElem).tempNodesGlobal=PHUTelem{indexPatch}(indexElem).nodesGlobal;
        PHUTelem{indexPatch}(indexElem).tempNodesGlobal(removeNodes)=0;
        allBasis=[allBasis,PHUTelem{indexPatch}(indexElem).tempNodesGlobal];
        
    end
end

%re-assign nodesGlobal
allBasis=nonzeros(allBasis);
allBasis=unique(allBasis);
numBasis=length(allBasis);
replacePattern=1:max(allBasis);
replacePattern(allBasis)=1:numBasis;
sizeBasis=numBasis;
for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        numBasis=length(PHUTelem{indexPatch}(indexElem).tempNodesGlobal);
        for ii=1:numBasis
            if PHUTelem{indexPatch}(indexElem).tempNodesGlobal(ii)~=0
                PHUTelem{indexPatch}(indexElem).tempNodesGlobal(ii)=replacePattern(PHUTelem{indexPatch}(indexElem).tempNodesGlobal(ii));
            end
        end
    end
end

% take type 2 basis from PHUTelemDE
allType2Basis=[];
for indexPatch=1:length(PHUTelemType2BasisInfo)
    for indexElem=1:length(PHUTelemType2BasisInfo{indexPatch})
        allType2Basis=[allType2Basis,PHUTelemType2BasisInfo{indexPatch}(indexElem).type2Nodes];
    end
end

%re-assign type 2 basis
allType2Basis=unique( allType2Basis);
numType2Basis=length(allType2Basis);
replacePattern=1:max(allType2Basis);
replacePattern(allType2Basis)=sizeBasis+1:sizeBasis+numType2Basis;
sizeBasis=sizeBasis+numType2Basis;
type2Basis=replacePattern(allType2Basis);
for indexPatch=1:length(PHUTelemDE)
    for indexElem=1:length(PHUTelemDE{indexPatch})
        PHUTelemType2BasisInfo{indexPatch}(indexElem).type2Nodes=replacePattern(PHUTelemType2BasisInfo{indexPatch}(indexElem).type2Nodes);
    end
end

% take additional basis from PHUTelemDE
allAdditionalBasis=[];
for indexPatch=1:length(PHUTelemDE)
    for indexElem=1:length(PHUTelemDE{indexPatch})
        additionalNodes=PHUTelemDE{indexPatch}(indexElem).additionalNodes;
        PHUTelemDE{indexPatch}(indexElem).additionalBasisC=PHUTelemDE{indexPatch}(indexElem).C(additionalNodes,:);
        PHUTelemDE{indexPatch}(indexElem).additionalBasis=PHUTelemDE{indexPatch}(indexElem).nodesGlobal(additionalNodes);
        allAdditionalBasis=[allAdditionalBasis,PHUTelemDE{indexPatch}(indexElem).additionalBasis];
    end
end

%re-assign additional basis
allAdditionalBasis=unique(allAdditionalBasis);
numAdditionalBasis=length(allAdditionalBasis);
replacePattern=1:max(allAdditionalBasis);
replacePattern(allAdditionalBasis)=sizeBasis+1:sizeBasis+numAdditionalBasis;
sizeBasis=sizeBasis+numAdditionalBasis;
additionalBasis=replacePattern(allAdditionalBasis);
for indexPatch=1:length(PHUTelemDE)
    for indexElem=1:length(PHUTelemDE{indexPatch})
        PHUTelemDE{indexPatch}(indexElem).additionalBasis=replacePattern(PHUTelemDE{indexPatch}(indexElem).additionalBasis);
    end
end

% assign type2 basis from PHUTelemType2BasisInfo to PHUTelem
for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        tempC=zeros(size(PHUTelem{indexPatch}(indexElem).modifiedC,1),(p+2)*(p+2));
        tempC(:,1:size(PHUTelem{indexPatch}(indexElem).modifiedC,2))=PHUTelem{indexPatch}(indexElem).modifiedC;
        PHUTelem{indexPatch}(indexElem).modifiedC=tempC;
        
        type2Locloc=find(PHUTelem{indexPatch}(indexElem).tempNodesGlobal==0);
        numType2Loc=length(type2Locloc);

        tempType2Basis=PHUTelemType2BasisInfo{indexPatch}(indexElem).type2Nodes;
        tempType2BasisC=PHUTelemType2BasisInfo{indexPatch}(indexElem).type2NodesC;
        numType2Basis=length(tempType2Basis);
        
        numBasis=length(PHUTelem{indexPatch}(indexElem).tempNodesGlobal);
        
        if numType2Basis<numType2Loc
            loc=type2Locloc(1:numType2Basis);
            PHUTelem{indexPatch}(indexElem).tempNodesGlobal(loc)= tempType2Basis;
            PHUTelem{indexPatch}(indexElem).modifiedC(loc,:)=tempType2BasisC;
            PHUTelem{indexPatch}(indexElem).polyDegree(loc)=p+1;
        else
            numAddNodes=numType2Basis-numType2Loc;
            loc=[type2Locloc,numBasis+1:numBasis+numAddNodes];
            PHUTelem{indexPatch}(indexElem).tempNodesGlobal(loc)= tempType2Basis;
            PHUTelem{indexPatch}(indexElem).modifiedC(loc,:)=tempType2BasisC;
            PHUTelem{indexPatch}(indexElem).polyDegree(loc)=p+1;
        end
    end
end

% assign additional basis from PHUTelemDE to PHUTelem
for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        
        additionalBasisloc=find(PHUTelem{indexPatch}(indexElem).tempNodesGlobal==0);
        numAdditionalLoc=length(additionalBasisloc);
        
        tempAdditionalBasis=PHUTelemDE{indexPatch}(indexElem).additionalBasis;
        tempAdditionalBasisC=PHUTelemDE{indexPatch}(indexElem).additionalBasisC;
        numAdditionalBasis=length(tempAdditionalBasis);
        
        numBasis=length(PHUTelem{indexPatch}(indexElem).tempNodesGlobal);
        
        if ~isempty(tempAdditionalBasis)
            if numAdditionalBasis<numAdditionalLoc
                
                loc=additionalBasisloc(1:numAdditionalBasis);
                PHUTelem{indexPatch}(indexElem).tempNodesGlobal(loc)= tempAdditionalBasis;
                PHUTelem{indexPatch}(indexElem).modifiedC(loc,:)=tempAdditionalBasisC;
                PHUTelem{indexPatch}(indexElem).polyDegree(loc)=p+1;
                
            else
                
                numAddNodes=numAdditionalBasis-numAdditionalLoc;
                loc=[additionalBasisloc,numBasis+1:numBasis+numAddNodes];
                PHUTelem{indexPatch}(indexElem).tempNodesGlobal(loc)= tempAdditionalBasis;
                PHUTelem{indexPatch}(indexElem).modifiedC(loc,:)=tempAdditionalBasisC;
                PHUTelem{indexPatch}(indexElem).polyDegree(loc)=p+1;
                
            end
        end
        
        
    end
end

% remove extraNodes & reset nodesGlobal
for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        extraLoc=find(PHUTelem{indexPatch}(indexElem).tempNodesGlobal==0);
        if ~isempty(extraLoc)
            PHUTelem{indexPatch}(indexElem).tempNodesGlobal(extraLoc)=[];
            PHUTelem{indexPatch}(indexElem).modifiedC(extraLoc,:)=[];
            PHUTelem{indexPatch}(indexElem).polyDegree(extraLoc)=[];
        end
        PHUTelem{indexPatch}(indexElem).nodesGlobal= PHUTelem{indexPatch}(indexElem).tempNodesGlobal;
    end
end

% figure
% plotPHTMesh_nodesGlobal3( PHUTelem,GIFTmesh,type2Basis,additionalBasis,p)
% pause

end

