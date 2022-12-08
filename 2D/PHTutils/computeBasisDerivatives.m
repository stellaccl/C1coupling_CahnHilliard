function [c_dRdu,c_dRdv] = computeBasisDerivatives(PHUTelem,GIFTmesh,indexPatch,indexElem,uref,vref,p,q )
%disp('in compute Basis derivatives')

[B_u,dB_u] = bernstein_basis(uref,p);
[B_v,dB_v] = bernstein_basis(vref,q);

R = zeros(1,1, (p+1)*(q+1));
dRdu = zeros(1,1, (p+1)*(q+1));
dRdv = zeros(1,1, (p+1)*(q+1));

basisCounter =0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        R(:,:,basisCounter) = B_u(:,i)*B_v(:,j)';
        dRdu(:,:,basisCounter) = dB_u(:,i)*B_v(:,j)';
        dRdv(:,:,basisCounter) = B_u(:,i)*dB_v(:,j)';
    end
end

numBasis=length(PHUTelem{indexPatch}(indexElem).nodesGlobal);

c_dRdu=zeros(1,length(numBasis));
c_dRdv=zeros(1,length(numBasis));

for index=1:numBasis
    
    c_dRdu(index) = PHUTelem{indexPatch}(indexElem).C(index,:)*dRdu(:);
    c_dRdv(index) = PHUTelem{indexPatch}(indexElem).C(index,:)*dRdv(:);
end

end

