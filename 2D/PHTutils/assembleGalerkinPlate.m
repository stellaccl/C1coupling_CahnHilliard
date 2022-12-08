function [ stiff, rhs ] = assembleGalerkinPlate( PHTelem, GIFTmesh, sizeBasis, p, q, Cmat, q0 )
%assembles the stiffness matrix and rhs (Galerkin method) for the plate
%problem
%uses GIFT mapping
%supports multipatches

%Gauss points
ngauss_x = p+2;
ngauss_y = q+2;
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


%initialize LHS stiffness matrix and RHS vector
stiff = sparse(sizeBasis,sizeBasis);
rhs = zeros(sizeBasis,1);

%assemble the stiffness matrix and RHS
domain_area = 0;
for patchIndex = 1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            xmin = PHTelem{patchIndex}(i).vertex(1);
            xmax = PHTelem{patchIndex}(i).vertex(3);
            ymin = PHTelem{patchIndex}(i).vertex(2);
            ymax = PHTelem{patchIndex}(i).vertex(4);
            
            %the jacobian of the transformation from [-1,1]x[-1,1] to
            %[xmin, xmax]x [ymin, ymax]
            scalefac = (xmax - xmin)*(ymax - ymin)/4;
            nument = size(PHTelem{patchIndex}(i).modifiedC,1);
            scrtx = PHTelem{patchIndex}(i).nodesGlobal(1:nument);
            %dscrtx = reshape([2*scrtx-1; 2*scrtx],1,2*nument);
            
            localstiff = zeros(nument, nument); %local stiffness
            localrhs = zeros(nument, 1);
            
            %loop over the ngauss_x x ngauss_y gauss points on each element
            for jj=1:ngauss_y
                for ii=1:ngauss_x
                    %evaluate the derivatives of the mapping from parameter
                    %space to physical space
                    
                    [ coord, dxdxi, d2xdxi2, dxdxi2] = paramMapPlate( GIFTmesh{patchIndex}, gauss_coord_x(ii), gauss_coord_y(jj), xmin, ymin, xmax, ymax);
%                     
%                     dxdxi
%                     d2xdxi2
%                     dxdxi2
%                     
% 
%                     
%                     plot(coord(1),coord(2),'+r')
%                     hold on
%                     drawnow
%                     
%                     pause
                    
                    R = (PHTelem{patchIndex}(i).C)*squeeze(Buv(ii,jj,:));
                    dRdx = (PHTelem{patchIndex}(i).C)*squeeze(dBdu(ii,jj,:));
                    dRdy = (PHTelem{patchIndex}(i).C)*squeeze(dBdv(ii,jj,:));
                    d2Rdx = (PHTelem{patchIndex}(i).C)*squeeze(d2Bdu(ii,jj,:));
                    d2Rdy = (PHTelem{patchIndex}(i).C)*squeeze(d2Bdv(ii,jj,:));
                    d2Rdxdy = (PHTelem{patchIndex}(i).C)*squeeze(d2Bdudv(ii,jj,:));
                    
                    %multiply by the jacobian of the transformation from reference
                    %space to the parameter space
                    dRdx = dRdx*2/(xmax-xmin);
                    dRdy = dRdy*2/(ymax-ymin);
                    d2Rdx = d2Rdx*(2/(xmax-xmin))^2;
                    d2Rdy = d2Rdy*(2/(ymax-ymin))^2;
                    d2Rdxdy = d2Rdxdy*(2/(xmax-xmin))*(2/(ymax-ymin));
                    
                    % Solve for first derivatives in global coordinates
                    dR = dxdxi\[dRdx';dRdy'];
                    d2R =  dxdxi2\([d2Rdx';d2Rdxdy';d2Rdy']-d2xdxi2*dR);
                    
%                     d2R
%                     pause
                    
                    B = zeros(3,nument);
                    J = det(dxdxi);
                    
                    B(1,:) = d2R(1,:);
                    B(2,:) = d2R(3,:);
                    B(3,:) = 2*d2R(2,:);
                               
                    localstiff = localstiff + B' * Cmat * B * scalefac * gauss_weight_x(ii).*gauss_weight_y(jj).*J;%
                    localrhs = localrhs + q0*R*scalefac*gauss_weight_x(ii).*gauss_weight_y(jj).*J;
                    domain_area = domain_area + scalefac * gauss_weight_x(ii).*gauss_weight_y(jj).*J;
                    
                end
            end
            stiff(scrtx, scrtx) = stiff(scrtx, scrtx) + localstiff;
            
%             scrtx
%             rhs(scrtx)
%             localrhs
            rhs(scrtx) = rhs(scrtx) + localrhs;
        end
    end
end


domain_area
domain_area_error = domain_area - 10^4