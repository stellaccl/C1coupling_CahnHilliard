function [PHUTelem,c1SizeBasis,type2Basis] = assignC1basis_degElev(PHUTelem,PHUTelemDE,c1SizeBasis,coefSol,solIndexCount,p,q)

numType2Basis=size(coefSol,2);
type2Basis=c1SizeBasis+1:c1SizeBasis+numType2Basis;
c1SizeBasis=c1SizeBasis+numType2Basis;

for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(indexElem).children)
            
            type2Loc=find(PHUTelem{indexPatch}(indexElem).remainNodes==0);
            
            if~isempty(type2Loc)
                tempNodes=PHUTelem{indexPatch}(indexElem).nodesGlobal;
                tempNodes=[tempNodes,type2Basis];
                PHUTelem{indexPatch}(indexElem).nodesGlobal=tempNodes;
                PHUTelem{indexPatch}(indexElem).polyDegree(end+1:end+numType2Basis)=p+1;
            end
            
        end
    end
end

[ PHUTelem ] = assignC1Function_degElev(PHUTelem,PHUTelemDE,coefSol,type2Basis,solIndexCount,p,q);

end

