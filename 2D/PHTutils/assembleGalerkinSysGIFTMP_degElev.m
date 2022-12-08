function [ stiff, rhs ] = assembleGalerkinSysGIFTMP_degElev( PHTelem, GIFTmesh, sizeBasis, p, q, Cmat )
%assembles the stiffness matrix and rhs (Galerkin method)
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
[B_u, dB_u] = bernstein_basis(gauss_coord_x,p);
[B_v, dB_v] = bernstein_basis(gauss_coord_y,q);

%1D bernstein polynomials evaluated at the Gauss points on the master
%element for the degree elevated basis
[B_uDE, dB_uDE] = bernstein_basis(gauss_coord_x,p+1);
[B_vDE, dB_vDE] = bernstein_basis(gauss_coord_y,q+1);


dBdu = zeros(ngauss_x, ngauss_y, (p+1)*(q+1));
dBdv = zeros(ngauss_x, ngauss_y, (p+1)*(q+1));

dBduDE = zeros(ngauss_x, ngauss_y, (p+2)*(q+2));
dBdvDE = zeros(ngauss_x, ngauss_y, (p+2)*(q+2));

%the derivatives of the 2D Bernstein polynomials at Gauss points on the
%master element
basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        dBdu(:,:,basisCounter) = dB_u(:,i)*B_v(:,j)';
        dBdv(:,:,basisCounter) = B_u(:,i)*dB_v(:,j)';
    end
end

%the derivatives of the 2D Bernstein degree elevated polynomials at Gauss points on the
%master element
basisCounter = 0;
for j=1:q+2
    for i=1:p+2
        basisCounter = basisCounter + 1;
        dBduDE(:,:,basisCounter) = dB_uDE(:,i)*B_vDE(:,j)';
        dBdvDE(:,:,basisCounter) = B_uDE(:,i)*dB_vDE(:,j)';
    end
end


%initialize LHS stiffness matrix and RHS vector
dim = 2; %the dimension of physical space
stiff = sparse(dim*sizeBasis,dim*sizeBasis);
rhs = zeros(dim*sizeBasis,1);

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
            dscrtx = reshape([2*scrtx-1; 2*scrtx],1,2*nument);
            
            localstiff = zeros(2*nument, 2*nument); %local stiffness
            
            %loop over the ngauss_x x ngauss_y gauss points on each element
            for jj=1:ngauss_y
                for ii=1:ngauss_x
                    %evaluate the derivatives of the mapping from parameter
                    %space to physical space
                    [coord, dxdxi] = paramMap( GIFTmesh{patchIndex}, gauss_coord_x(ii), gauss_coord_y(jj), xmin, ymin, xmax, ymax);
                    J = det(dxdxi);
                    B = zeros(2*nument,3);

%                     plot(coord(1),coord(2),'+r')
%                     hold on
%                     drawnow
                    
                    for t=1:nument
                        if PHTelem{patchIndex}(i).polyDegree(t)==p
                            dRdx = (PHTelem{patchIndex}(i).modifiedC(t,1:(p+1)*(q+1)))*squeeze(dBdu(ii,jj,:));
                            dRdy = (PHTelem{patchIndex}(i).modifiedC(t,1:(p+1)*(q+1)))*squeeze(dBdv(ii,jj,:));
                        else
                          
                            dRdx = (PHTelem{patchIndex}(i).modifiedC(t,:))*squeeze(dBduDE(ii,jj,:));
                            dRdy = (PHTelem{patchIndex}(i).modifiedC(t,:))*squeeze(dBdvDE(ii,jj,:));
                        end
                        
                        %multiply by the jacobian of the transformation from reference
                        %space to the parameter space
                        dRdx = dRdx*2/(xmax-xmin);
                        dRdy = dRdy*2/(ymax-ymin);
                        
                        % Solve for first derivatives in global coordinates
                        dR = dxdxi\[dRdx';dRdy'];
                                               
                        B(2*t-1,1) = dR(1);
                        B(2*t,2) = dR(2);
                        B(2*t-1,3) = dR(2);
                        B(2*t,3) = dR(1);
                    end
%                     B
%                     pause
                    %TODO: implement non-zero volume force
                    localstiff = localstiff + B * Cmat * B' * scalefac * gauss_weight_x(ii).*gauss_weight_y(jj).*J;
                    domain_area = domain_area + scalefac * gauss_weight_x(ii).*gauss_weight_y(jj).*J;
                    
                end
            end
            stiff(dscrtx, dscrtx) = stiff(dscrtx, dscrtx) + localstiff;
        end
    end
end


domain_area
domain_area_error = domain_area - (16-pi/4)
