function [PHUTelem,c1SizeBasis,allType2Basis] = zipConformingC1_localization(PHUTelem,c1SizeBasis,allType2Basis,coefSol,solIndexCount,p,q )

numType2Basis=size(coefSol,2);
type2Basis=c1SizeBasis+1:c1SizeBasis+numType2Basis;
allType2Basis=[allType2Basis,type2Basis];
c1SizeBasis=c1SizeBasis+numType2Basis;
for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(indexElem).children)
            
            type2Loc=find(PHUTelem{indexPatch}(indexElem).solIndex~=0);
            
            if~isempty(type2Loc)
                
                tempNodes=PHUTelem{indexPatch}(indexElem).C1NodesGlobal;
                tempNodes=[tempNodes,type2Basis];
                PHUTelem{indexPatch}(indexElem).C1NodesGlobal=tempNodes;
                PHUTelem{indexPatch}(indexElem).polyDegree(end+1:end+numType2Basis)=p;
            end
            
        end
    end
end

[ PHUTelem ] = assignC1Function(PHUTelem,coefSol,type2Basis,solIndexCount,p,q);
end

