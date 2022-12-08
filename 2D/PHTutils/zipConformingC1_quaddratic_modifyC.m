function[PHUTelem,sizeBasis,type2Basis] = zipConformingC1_quaddratic_modifyC(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,coefSol,solIndexCount,p,q)


numType2Basis=size(coefSol,2);

type2basisNodes=sizeBasis+1:sizeBasis+numType2Basis;

[PHUTelem,type2Basis,sizeBasis] = assignNewNodesGlobal_quadratic(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,type2basisNodes,p,q);


[PHUTelem] = modifyC_quadratic(PHUTelem,coefSol,numType2Basis,type2Basis,solIndexCount,p,q);

end

