function [PHUTelem,PHUTelemDE,c1SizeBasis ] = assignAdditionalNodes_degElev( PHUTelem,PHUTelemDE,c1SizeBasis,p,q)
%take additional basis from PHUTelemDE
%reassig nodes global (aal equal zeros except for additional basis)

allAdditionalBasis=[];
for indexPatch=1:length(PHUTelemDE)
    for indexElem=1:length(PHUTelemDE{indexPatch})

         PHUTelemDE{indexPatch}(indexElem).remainNodes=PHUTelemDE{indexPatch}(indexElem).nodesGlobal;
        
        additionalNodes=PHUTelemDE{indexPatch}(indexElem).additionalNodes;
        PHUTelemDE{indexPatch}(indexElem).additionalBasisC=PHUTelemDE{indexPatch}(indexElem).C(additionalNodes,:);
        PHUTelemDE{indexPatch}(indexElem).additionalBasis=PHUTelemDE{indexPatch}(indexElem).nodesGlobal(additionalNodes);
        
        allNodes=1:length(PHUTelemDE{indexPatch}(indexElem).nodesGlobal);
        allAdditionalBasis=[allAdditionalBasis,PHUTelemDE{indexPatch}(indexElem).additionalBasis];
        
        zerosNodes= setdiff(allNodes,additionalNodes);
        PHUTelemDE{indexPatch}(indexElem).remainNodes(zerosNodes)=0;
        
    end
end

allAdditionalBasis=unique(allAdditionalBasis);
numAdditionalBasis=length(allAdditionalBasis);
replacePattern=1:max(allAdditionalBasis);
replacePattern(allAdditionalBasis)=c1SizeBasis+1:c1SizeBasis+numAdditionalBasis;

allAdditionalBasis=replacePattern(allAdditionalBasis);

c1SizeBasis=c1SizeBasis+numAdditionalBasis;
for indexPatch=1:length(PHUTelemDE)
    for indexElem=1:length(PHUTelemDE{indexPatch})
        PHUTelemDE{indexPatch}(indexElem).additionalBasis=replacePattern(PHUTelemDE{indexPatch}(indexElem).additionalBasis);
        
        for i=1:length(PHUTelemDE{indexPatch}(indexElem).remainNodes)
            if PHUTelemDE{indexPatch}(indexElem).remainNodes(i)~=0
                PHUTelemDE{indexPatch}(indexElem).remainNodes(i)=replacePattern(PHUTelemDE{indexPatch}(indexElem).remainNodes(i));
            end
        end
    end
end

%assign into PHTelem the additional nodes & corresponding modifiedC
for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        
        tempAdditionalBasis=PHUTelemDE{indexPatch}(indexElem).additionalBasis;
        tempAdditionalBasisC=PHUTelemDE{indexPatch}(indexElem).additionalBasisC;
        numAdditionalBasis=length(tempAdditionalBasis);
        if ~isempty(tempAdditionalBasis)
            tempNodes=PHUTelem{indexPatch}(indexElem).nodesGlobal;
            tempNodes=[tempNodes,tempAdditionalBasis];
            PHUTelem{indexPatch}(indexElem).nodesGlobal=tempNodes;
            
            PHUTelem{indexPatch}(indexElem).modifiedC(end+1:end+numAdditionalBasis,:)=tempAdditionalBasisC;
            PHUTelem{indexPatch}(indexElem).polyDegree(end+1:end+numAdditionalBasis)=p+1;  
             
        end
        
        
    end
end

end

