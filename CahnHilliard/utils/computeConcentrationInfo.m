function  [localC, dLocalC, del2_c,localC_t]=computeConcentrationInfo(PHTelem, GIFTmesh,  p, q, C,Cdot,uref,vref,patchIndex,elemIndex)
%compute concentration C , its derivatives, and C_t at a given quadrature
%point (uref vref)

numPts = 1; %number of plot points to use on each edge

%1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u,dB_u,ddB_u] = bernstein_basis(uref,p);
[B_v, dB_v,ddB_v] = bernstein_basis(vref,q);

Buv = zeros(numPts, numPts, (p+1)*(q+1));
dBdu = zeros(numPts, numPts, (p+1)*(q+1));
dBdv = zeros(numPts, numPts, (p+1)*(q+1));
d2Bdu = zeros(numPts, numPts, (p+1)*(q+1));
d2Bdv = zeros(numPts, numPts, (p+1)*(q+1));
d2Bdudv = zeros(numPts, numPts, (p+1)*(q+1));
%the derivatives of the 2D Bernstein polynomials at Gauss points on the
%master element
basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        Buv(:,:,basisCounter) = B_u(:,i)*B_v(:,j)';
        dBdu(:,:,basisCounter) = dB_u(:,i)*B_v(:,j)';
        dBdv(:,:,basisCounter) = B_u(:,i)*dB_v(:,j)';
        d2Bdu(:,:,basisCounter) = ddB_u(:,i)*B_v(:,j)';
        d2Bdv(:,:,basisCounter) = B_u(:,i)*ddB_v(:,j)';
        d2Bdudv(:,:,basisCounter) = dB_u(:,i)*dB_v(:,j)';
    end
end


xmin = PHTelem{patchIndex}(elemIndex).vertex(1);
xmax = PHTelem{patchIndex}(elemIndex).vertex(3);
ymin = PHTelem{patchIndex}(elemIndex).vertex(2);
ymax = PHTelem{patchIndex}(elemIndex).vertex(4);

%the jacobian of the transformation from [-1,1]x[-1,1] to
%[xmin, xmax]x [ymin, ymax]
scalefac = (xmax - xmin)*(ymax - ymin)/4;

% globalNodes=PHTelem{patchIndex}(i).nodesGlobal;
tempC=C(PHTelem{patchIndex}(elemIndex).nodes);
tempCdot=Cdot(PHTelem{patchIndex}(elemIndex).nodes);

R = (PHTelem{patchIndex}(elemIndex).C)*squeeze(Buv(1,1,:));
dRdx = (PHTelem{patchIndex}(elemIndex).C)*squeeze(dBdu(1,1,:));
dRdy = (PHTelem{patchIndex}(elemIndex).C)*squeeze(dBdv(1,1,:));
d2Rdx = (PHTelem{patchIndex}(elemIndex).C)*squeeze(d2Bdu(1,1,:));
d2Rdy = (PHTelem{patchIndex}(elemIndex).C)*squeeze(d2Bdv(1,1,:));
d2Rdxdy = (PHTelem{patchIndex}(elemIndex).C)*squeeze(d2Bdudv(1,1,:));

%multiply by the jacobian of the transformation from reference
%space to the parameter space
dRdx = dRdx*2/(xmax-xmin);
dRdy = dRdy*2/(ymax-ymin);
d2Rdx = d2Rdx*(2/(xmax-xmin))^2;
d2Rdy = d2Rdy*(2/(ymax-ymin))^2;
d2Rdxdy = d2Rdxdy*(2/(xmax-xmin))*(2/(ymax-ymin));

[~, dxdxi, d2xdxi2, dxdxi2] = paramMapPlate( GIFTmesh{patchIndex},uref,vref, xmin, ymin, xmax, ymax);
%J = det(dxdxi);
dR = dxdxi\[dRdx'; dRdy'];
ddR=[d2Rdx';d2Rdxdy';d2Rdy'];
ddR = dxdxi2\(ddR - d2xdxi2*dR);

localC=R'*tempC;
dcdu=dR(1,:)*tempC;
dcdv=dR(2,:)*tempC;
d2cdu=ddR(1,:)*tempC;
d2cdv=ddR(3,:)*tempC;
localC_t=R'*tempCdot;
del2_c = d2cdu+d2cdv; %laplacian
dLocalC=[dcdu,dcdv];



end