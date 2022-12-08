function [PHUTelem] = modifyC_quadratic(PHUTelem,coefSol,numType2Basis,type2Basis,solIndexCount,p,q)
%compute modifiedC = sol* associated bazier coef (method 2)
%assign modifiedC to relevent basis

%for each sol Index, assign the associanted basis (patch, elem, nodes of the basis)
for indexPatch=1:length(PHUTelem)
    
    for indexElem=1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(indexElem).children)
    
            PHUTelem{indexPatch}(indexElem).modifiedC=PHUTelem{indexPatch}(indexElem).C;
            
            type2Loc=find(PHUTelem{indexPatch}(indexElem).solIndex~=0);
            
            if ~isempty(type2Loc)
                
                solIndexInfo=zeros(solIndexCount,3);
                
                solIndex=PHUTelem{indexPatch}(indexElem).solIndex;
                for indexBasis=1:length(solIndex)
                    if solIndex(indexBasis)~=0
                        solIndexInfo(solIndex(indexBasis),:)=[indexPatch,indexElem,indexBasis];
                    end
                end
                
                mBezierCoef=zeros((p+1)*(q+1),solIndexCount);
                for indexColum=1:solIndexCount
                    if all(solIndexInfo(indexColum,:))
                        patch=solIndexInfo(indexColum,1);
                        elem=solIndexInfo(indexColum,2);
                        basis=solIndexInfo(indexColum,3);
                        tempCoef=PHUTelem{patch}(elem).C(basis,:);
                        mBezierCoef(:,indexColum)=tempCoef;
                    end
                end
                
                mModifiedC=zeros((p+1)*(q+1),numType2Basis);
                for indexBasis=1:numType2Basis
                    sol=coefSol(:,indexBasis);
                    mModifiedC(:,indexBasis)=mBezierCoef*sol;
                end
                
                [~,loc]=ismember(type2Basis,PHUTelem{indexPatch}(indexElem).nodesGlobal);
                type2NodesIndex=loc;
                
                
                if ~isempty(PHUTelem{indexPatch}(indexElem).extraNodes)
                    for ii=1:numType2Basis
                        localIndex =  type2NodesIndex(ii);
                        PHUTelem{indexPatch}(indexElem).modifiedC(localIndex,:) = mModifiedC(:,ii)';
                    end
                    extraNodes=PHUTelem{indexPatch}(indexElem).extraNodes;
                    PHUTelem{indexPatch}(indexElem).modifiedC(extraNodes,:)=[];
                else
                    for ii=1:numType2Basis
                        localIndex =  type2NodesIndex(ii);
                        PHUTelem{indexPatch}(indexElem).modifiedC(localIndex,1:size(mModifiedC(:,ii)',2)) = mModifiedC(:,ii)';
                    end
                end
            end
            
        end
    end
end


end




