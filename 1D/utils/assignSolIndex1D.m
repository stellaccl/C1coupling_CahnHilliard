function [PHTelem,solIndexCount] = assignSolIndex1D( PHTelem,p )
%only for 2 patch domain, coupling the right side of patch 1 with left side
%of patch 2

rightNodes=[p,p+1];
leftNodes=[1,2];

%select basis that introduce C0 continuity
indexPatch=1;
PHTelem{indexPatch}(end).nodesGlobal
tempNodes1=PHTelem{indexPatch}(end).nodesGlobal(rightNodes);

indexPatch=2;
tempNodes2=PHTelem{indexPatch}(1).nodesGlobal(leftNodes);
selectedBasis=unique([tempNodes1,tempNodes2]);

solIndex=1:length(selectedBasis);
solIndexCount=length(selectedBasis);

replacePatter=zeros(1,max(PHTelem{2}(end).nodesGlobal));
replacePatter(selectedBasis)=solIndex;

for indexPatch=1:length(PHTelem)
    for indexElem=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(indexElem).children)
            PHTelem{indexPatch}(indexElem).solIndex=replacePatter(PHTelem{indexPatch}(indexElem).nodesGlobal);
        end
    end
end


end

