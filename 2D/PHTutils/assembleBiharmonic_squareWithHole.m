function [stiff, rhs ] = assembleBiharmonic_squareWithHole(PHTelem,GIFTmesh, sizeBasis, p, q )
%use nodesGlobal
%use modifiedC
%vectorise stiffness matrix
%Gauss points
ngauss_x = p+1;
ngauss_y = q+1;
[gauss_weight_x, gauss_coord_x] = quadrature( ngauss_x, 'GAUSS', 1 );
[gauss_weight_y, gauss_coord_y] = quadrature( ngauss_y, 'GAUSS', 1 );

%take the transpose so that they are in the format expected by
%bernstein_basis
gauss_coord_x = gauss_coord_x';
gauss_coord_y = gauss_coord_y';

%1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u, dB_u,ddB_u] = bernstein_basis(gauss_coord_x,p);
[B_v, dB_v,ddB_v] = bernstein_basis(gauss_coord_y,q);

Buv = zeros(ngauss_x, ngauss_y, (p+1)*(q+1));
dBdu = zeros(ngauss_x, ngauss_y, (p+1)*(q+1));
dBdv = zeros(ngauss_x, ngauss_y, (p+1)*(q+1));
d2Bdu = zeros(ngauss_x, ngauss_y, (p+1)*(q+1));
d2Bdv = zeros(ngauss_x, ngauss_y, (p+1)*(q+1));
d2Bdudv = zeros(ngauss_x, ngauss_y, (p+1)*(q+1));

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



stiff = sparse(sizeBasis,sizeBasis);
rhs= zeros(sizeBasis,1);

for patchIndex = 1:length(PHTelem)
    for elemIndex=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(elemIndex).children)
            
            xmin = PHTelem{patchIndex}(elemIndex).vertex(1);
            xmax = PHTelem{patchIndex}(elemIndex).vertex(3);
            ymin = PHTelem{patchIndex}(elemIndex).vertex(2);
            ymax = PHTelem{patchIndex}(elemIndex).vertex(4);
            
            scalefac = (xmax - xmin)*(ymax - ymin)/4;
            nument = size(PHTelem{patchIndex}(elemIndex).modifiedC,1);
            scrtx = PHTelem{patchIndex}(elemIndex).nodesGlobal(1:nument);
            
            localStiff = zeros(nument, nument);
            localRhs = zeros(nument, 1);
            
            for jj=1:ngauss_y
                for ii=1:ngauss_x
                    
               
                    R = (PHTelem{patchIndex}(elemIndex).modifiedC)*squeeze(Buv(ii,jj,:));
                    dRdx = (PHTelem{patchIndex}(elemIndex).modifiedC)*squeeze(dBdu(ii,jj,:));
                    dRdy = (PHTelem{patchIndex}(elemIndex).modifiedC)*squeeze(dBdv(ii,jj,:));
                    
                    d2Rdx = (PHTelem{patchIndex}(elemIndex).modifiedC)*squeeze(d2Bdu(ii,jj,:));
                    d2Rdy = (PHTelem{patchIndex}(elemIndex).modifiedC)*squeeze(d2Bdv(ii,jj,:));
                    d2Rdxdy = (PHTelem{patchIndex}(elemIndex).modifiedC)*squeeze(d2Bdudv(ii,jj,:));
                                        
                    %multiply by the jacobian of the transformation from reference
                    %space to the parameter space
                    dRdx = dRdx*2/(xmax-xmin);
                    dRdy = dRdy*2/(ymax-ymin);
                    
                    d2Rdx = d2Rdx*(2/(xmax-xmin))^2;
                    d2Rdy = d2Rdy*(2/(ymax-ymin))^2;
                    
                    d2Rdxdy = d2Rdxdy*(2/(xmax-xmin))*(2/(ymax-ymin));

                    [coords, dxdxi, d2xdxi2, dxdxi2] = paramMapPlate( GIFTmesh{patchIndex},gauss_coord_x(ii), gauss_coord_y(jj), xmin, ymin, xmax, ymax);
                   
                    dR = dxdxi\[dRdx'; dRdy'];
                    ddR=[d2Rdx';d2Rdxdy';d2Rdy'];
                    ddR = dxdxi2\(ddR - d2xdxi2*dR);
                    
                    J = abs(det(dxdxi));
                    
                    ddRx= ddR(1,:);
                    ddRy= ddR(3,:);

                    part1=ddRx'*(ddRx+ddRy);
                    part2=ddRy'*(ddRx+ddRy);
                    
                    tempValue=part1+part2;
                    localStiff=localStiff+tempValue*scalefac*gauss_weight_x(ii).*gauss_weight_y(jj).*J;
                    
                    [~,~,d4udx4,d4udy4,d4udx2dy2]=computeDerivativesSquareWithHoleBiharmonic(coords(1),coords(2));
                
                    f=d4udx4+d4udy4+2*d4udx2dy2;
                    
                    localRhsTemp=R*f;
                    
                    localRhsTemp=localRhsTemp.*scalefac.*gauss_weight_x(ii).*gauss_weight_y(jj).*J;
                    
                    localRhs = localRhs + localRhsTemp;
                end
            end
            
            stiff(scrtx,scrtx) = stiff(scrtx,scrtx) + localStiff;
            rhs(scrtx) =  rhs(scrtx) + localRhs;
        end
    end
end


end

