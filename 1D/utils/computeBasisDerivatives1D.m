function [c_dRdu] =computeBasisDerivatives1D( PHUTelem,GIFTmesh,indexPatch,indexElem,uref,p)

[B_u,dB_u] = bernstein_basis(uref,p);
numBasis=length(PHUTelem{indexPatch}(indexElem).nodesGlobal);
c_dRdu=zeros(1,length(numBasis));

xmin = PHUTelem{indexPatch}(indexElem).vertex(1);
xmax = PHUTelem{indexPatch}(indexElem).vertex(2);
[~,dxdxi] = paramMap1D(GIFTmesh{indexPatch}, uref, xmin, xmax );

for index=1:numBasis
    tempDRdu = PHUTelem{indexPatch}(indexElem).C(index,:)*dB_u(:);
%    c_dRdu(index) = PHUTelem{indexPatch}(indexElem).C(index,:)*dB_u(:);
    c_dRdu(index) = dxdxi\tempDRdu;
end

end

