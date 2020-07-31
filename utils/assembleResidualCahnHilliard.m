function [ residual ] = assembleResidualCahnHilliard( PHTelem, GIFTmesh, sizeBasis, p, q,C,Cdot,theta,alpha)
%assembles residual (CahnHilliard equation)
%disp('in assemble Residual')
% C
% pause
% Cdot
% pause
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


%initialize residual vector
residual= zeros(sizeBasis,1);

for patchIndex = 1:length(PHTelem)
    for elemIndex=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(elemIndex).children)
            xmin = PHTelem{patchIndex}(elemIndex).vertex(1);
            xmax = PHTelem{patchIndex}(elemIndex).vertex(3);
            ymin = PHTelem{patchIndex}(elemIndex).vertex(2);
            ymax = PHTelem{patchIndex}(elemIndex).vertex(4);
            
            %the jacobian of the transformation from [-1,1]x[-1,1] to
            %[xmin, xmax]x [ymin, ymax]
            scalefac = (xmax - xmin)*(ymax - ymin)/4;
            nument = size(PHTelem{patchIndex}(elemIndex).C,1);
            scrtx = PHTelem{patchIndex}(elemIndex).nodes(1:nument);
            %dscrtx = reshape([2*scrtx-1; 2*scrtx],1,2*nument);
            
            localResidual = zeros(nument, 1);
%             clc
%             elemIndex
            %loop over the ngauss_x x ngauss_y gauss points on each element
            for jj=1:ngauss_y
                for ii=1:ngauss_x
                    %evaluate the derivatives of the mapping from parameter
                    %space to physical space
                    
                    R = (PHTelem{patchIndex}(elemIndex).C)*squeeze(Buv(ii,jj,:));
                    dRdx = (PHTelem{patchIndex}(elemIndex).C)*squeeze(dBdu(ii,jj,:));
                    dRdy = (PHTelem{patchIndex}(elemIndex).C)*squeeze(dBdv(ii,jj,:));
                    d2Rdx = (PHTelem{patchIndex}(elemIndex).C)*squeeze(d2Bdu(ii,jj,:));
                    d2Rdy = (PHTelem{patchIndex}(elemIndex).C)*squeeze(d2Bdv(ii,jj,:));
                    d2Rdxdy = (PHTelem{patchIndex}(elemIndex).C)*squeeze(d2Bdudv(ii,jj,:));
                    
                    %multiply by the jacobian of the transformation from reference
                    %space to the parameter space
                    dRdx = dRdx*2/(xmax-xmin);
                    dRdy = dRdy*2/(ymax-ymin);
                    d2Rdx = d2Rdx*(2/(xmax-xmin))^2;
                    d2Rdy = d2Rdy*(2/(ymax-ymin))^2;
                    d2Rdxdy = d2Rdxdy*(2/(xmax-xmin))*(2/(ymax-ymin));
                    
                    [coord, dxdxi, d2xdxi2, dxdxi2] = paramMapPlate( GIFTmesh{patchIndex},gauss_coord_x(ii), gauss_coord_y(jj), xmin, ymin, xmax, ymax);
                   % coord
                    %                     ii
                    %                     jj
                    
                    J = det(dxdxi);
                    dR = dxdxi\[dRdx'; dRdy'];
                    ddR=[d2Rdx';d2Rdxdy';d2Rdy'];
                    ddR = dxdxi2\(ddR - d2xdxi2*dR);
                    
                    uref=gauss_coord_x(ii);
                    vref=gauss_coord_y(jj);
                    
                    [localC, dLocalC, del2_c,localC_t]=computeConcentrationInfo(PHTelem, GIFTmesh,  p, q, C,Cdot,uref,vref,patchIndex,elemIndex);
                    [mu,dmu,d2mu ] = chemicalPotential( localC, theta,alpha );
                    [M,dM,d2M] = mobility( localC);
                    %localC
                   
                    %d2M
%                     [del2_c;
%                         dLocalC(1);
%                         dLocalC(2)]
%                      M
%                     dM
%                     R
%                     dR
                    
%                     if elemIndex==5
%                         pause
%                     end
                    t1=M*dmu+dM*del2_c;
                    
                    localResidualTemp=R*localC_t+( dR(1,:)'*dLocalC(1)+dR(2,:)'*dLocalC(2))*t1+(ddR(1,:)'+ddR(3,:)')*M*del2_c;
                    %pause
                    localResidualTemp=localResidualTemp.*scalefac.*gauss_weight_x(ii).*gauss_weight_y(jj).*J;
                    %gauss_weight_x(ii).*gauss_weight_y(jj)
                    localResidual = localResidual + localResidualTemp;
                    %localResidual
                    
                end
            end
            %             localResidual
            %             if elemIndex==2
            %                 pause
            %             end
            residual(scrtx) =  residual(scrtx) + localResidual;
        end
    end
end
% residual
% pause
% 


end

