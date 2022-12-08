function [ R,dRdu,  dRdv ] = tensorProductBernsteinPolynomial(uref,vref,p,q )

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
end

