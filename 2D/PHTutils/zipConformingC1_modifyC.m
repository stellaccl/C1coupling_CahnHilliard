function [ PHUTelem,sizeBasis,type2Basis] = zipConformingC1_modifyC( PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,coefSol,p,q )
[boundaryInfo,~,~,~,~,patchEdgeInfo] = createBoundaryInfo3(PHUTelem, patchBoundaries );

numType2Basis=size(coefSol,2)

type2basisNodes=sizeBasis+1:sizeBasis+numType2Basis;
[PHUTelem ,type2Basis,sizeBasis] = assignNewNodesGlobal(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,type2basisNodes,p,q);
[PHUTelem] = modifyC(PHUTelem,coefSol,numType2Basis,patchEdgeInfo,boundaryInfo,type2Basis,p,q);


end

