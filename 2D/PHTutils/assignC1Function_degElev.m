function [ PHUTelem ] = assignC1Function_degElev(PHUTelem,PHUTelemDE,coefSol,type2Basis,solIndexCount,p,q)

numType2Basis=length(type2Basis);

%assign new basis
for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(indexElem).children)

            type2Loc=find(PHUTelemDE{indexPatch}(indexElem).solIndex~=0);
            
            if~isempty(type2Loc)
                
                solIndexInfo=zeros(solIndexCount,3);
                
                solIndex=PHUTelemDE{indexPatch}(indexElem).solIndex;
                for indexBasis=1:length(solIndex)
                    if solIndex(indexBasis)~=0
                        solIndexInfo(solIndex(indexBasis),:)=[indexPatch,indexElem,indexBasis];
                    end
                end
                
                mBezierCoef=zeros((p+2)*(q+2),solIndexCount);
                for indexColum=1:solIndexCount
                    if all(solIndexInfo(indexColum,:))
                        patch=solIndexInfo(indexColum,1);
                        elem=solIndexInfo(indexColum,2);
                        basis=solIndexInfo(indexColum,3);
                        tempCoef=PHUTelemDE{patch}(elem).C(basis,:);
                        mBezierCoef(:,indexColum)=tempCoef;
                    end
                end
                
                mModifiedC=zeros((p+2)*(q+2),numType2Basis);
                for indexBasis=1:numType2Basis
                    sol=coefSol(:,indexBasis);
                    mModifiedC(:,indexBasis)=mBezierCoef*sol;
                end
                
                [~,loc]=ismember(type2Basis,PHUTelem{indexPatch}(indexElem).nodesGlobal);
                type2NodesIndex=loc;
                
%                 indexPatch
%                 indexElem
%                 
%                 PHUTelem{indexPatch}(indexElem).nodesGlobal
%                 type2Basis
                
                
                
                for ii=1:numType2Basis
                    localIndex =  type2NodesIndex(ii);
                    
                    
                    PHUTelem{indexPatch}(indexElem).modifiedC(localIndex,1:size(mModifiedC(:,ii)',2)) = mModifiedC(:,ii)';
                end

                
            end
     
        end
    end
    
end

