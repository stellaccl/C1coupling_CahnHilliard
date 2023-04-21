function [ PHTelem,type2Basis,sizeBasis] = zipConformingC1_1D(PHTelem,GIFTmesh,p,sizeBasis,solIndexCount)
% only for two patch domain
% join right side of patch 1 with left side of patch 2
indexPatch=1;
indexElem=length(PHTelem{indexPatch});
uref=1;
[dRduPatch1] =computeBasisDerivatives1D( PHTelem,GIFTmesh,indexPatch,indexElem,uref,p);

indexPatch=2;
indexElem=1;
uref=-1;
[dRduPatch2] =computeBasisDerivatives1D( PHTelem,GIFTmesh,indexPatch,indexElem,uref,p);

%m=[dRduPatch1(2) dRduPatch1(3)-dRduPatch2(1) -dRduPatch2(2)];
m=[dRduPatch1(end-1) dRduPatch1(end)-dRduPatch2(1) -dRduPatch2(2)];
coefSol=null(m);

numType2Basis=size(coefSol,2);
type2Basis=sizeBasis+1:sizeBasis+numType2Basis;

[PHTelem ,type2Basis,sizeBasis] = assignNewNodesGlobal1D( PHTelem,type2Basis,p);

[PHTelem] = modifyC_1D(PHTelem,coefSol,solIndexCount ,p);


end

