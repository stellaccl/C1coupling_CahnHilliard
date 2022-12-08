function [ stiff, rhs] = assembleGalerkin_shellC1( PHTelem, GIFTmesh, sizeBasis, p, q, nu, memStiff, benStiff  )
%assembles the stiffness matrix and rhs (Galerkin method) for the
%pinchedCylinder problem
%uses GIFT mapping

%Gauss points
ngauss_x = p+3;
ngauss_y = q+3;
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
stiff = sparse(3*sizeBasis,3*sizeBasis);
rhs = zeros(3*sizeBasis,1);

%assemble the stiffness matrix and RHS

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
            sctr = PHTelem{patchIndex}(i).nodesGlobal(1:nument);
            
            
            %tsctrx = reshape([3*sctr-2, 3*sctr-1, 3*sctr],1,3*nument)
            tsctrx = zeros(1, 3*nument);
            tsctrx(1:3:3*nument) = 3*sctr-2;
            tsctrx(2:3:3*nument) = 3*sctr-1;
            tsctrx(3:3:3*nument) = 3*sctr;

            
            localstiff = zeros(3*nument, 3*nument); %local stiffness
            %localrhs = zeros(nument, 1);
            
            %loop over the ngauss_x x ngauss_y gauss points on each element
            for jj=1:ngauss_y
                for ii=1:ngauss_x
                    %evaluate the derivatives of the mapping from parameter
                    %space to physical space
                    %gauss_coord_x(ii)
                    %gauss_coord_y(jj)
                    [coord, dxdxi, d2xdxi2] = paramMapShell( GIFTmesh{patchIndex}, gauss_coord_x(ii), gauss_coord_y(jj), xmin, ymin, xmax, ymax);
                    %plot3(coord(:,1),coord(:,2), coord(:,3), '+r')
                    hold on
                    
                    R = (PHTelem{patchIndex}(i).modifiedC)*squeeze(Buv(ii,jj,:));
                    dRdx = (PHTelem{patchIndex}(i).modifiedC)*squeeze(dBdu(ii,jj,:));
                    dRdy = (PHTelem{patchIndex}(i).modifiedC)*squeeze(dBdv(ii,jj,:));
                    d2Rdx = (PHTelem{patchIndex}(i).modifiedC)*squeeze(d2Bdu(ii,jj,:));
                    d2Rdy = (PHTelem{patchIndex}(i).modifiedC)*squeeze(d2Bdv(ii,jj,:));
                    d2Rdxdy = (PHTelem{patchIndex}(i).modifiedC)*squeeze(d2Bdudv(ii,jj,:));
                  %  pause
                    
                    %multiply by the jacobian of the transformation from reference
                    %space to the parameter space
                    dRdx = dRdx*2/(xmax-xmin);
                    dRdy = dRdy*2/(ymax-ymin);
                    d2Rdx = d2Rdx*(2/(xmax-xmin))^2;
                    d2Rdy = d2Rdy*(2/(ymax-ymin))^2;
                    d2Rdxdy = d2Rdxdy*(2/(xmax-xmin))*(2/(ymax-ymin));
   
                    a1    = dxdxi(1,:);
                    a2    = dxdxi(2,:);
                    a3    = cross(a1,a2);
                    norma = norm(a3);
                    a3    = a3/norma;
                    J1    = norma;
                    
                    a11   = d2xdxi2(1,:);
                    a22   = d2xdxi2(3,:);
                    a12   = d2xdxi2(2,:);
                    
                    % dot products of ai and ei
                    a1e1  = a1(1);
                    a1e2  = a1(2);
                    a1e3  = a1(3);
                    a2e1  = a2(1);
                    a2e2  = a2(2);
                    a2e3  = a2(3);
                    
                    % R_I,2*a1 + R_I,1*a2 for all shape functions
                    noBasis = length(R);
                    dRIa    = zeros(3,noBasis);
                    for indexBasis=1:noBasis
                        dRIa(:,indexBasis) = dRdy(indexBasis)*a1 + dRdx(indexBasis)*a2;
                    end

                    % compute the constitutive matrix C
                    a_11 = dot(a1,a1);
                    a_12 = dot(a1,a2);
                    a_21 = dot(a2,a1);
                    a_22 = dot(a2,a2);
                    
                    aa1 = [a_11 a_21; a_12 a_22 ] \ [1;0];
                    aa2 = [a_11 a_21; a_12 a_22 ] \ [0;1];
                    
                    au11 = aa1(1);
                    au12 = aa1(2);
                    au22 = aa2(2);
                    
                    %                     Cmat = [au11^2 nu*au11*au22+(1-nu)*au12^2 au11*au12;
                    %                         au11*au12 au22*au12 0.5*((1-nu)*au11*au22+(1+nu)*au12^2);
                    %                         nu*au11*au22+(1-nu)*au12^2 au22^2 au22*au12];
                    
                    %ori
                    Cmat = [au11^2 nu*au11*au22+(1-nu)*au12^2 au11*au12;
                        nu*au11*au22+(1-nu)*au12^2 au22^2 au22*au12;
                        au11*au12 au22*au12 0.5*((1-nu)*au11*au22+(1+nu)*au12^2)];
         
                    Bmem = zeros(3,3*noBasis);
                    Bben = zeros(3,3*noBasis);
                    
                    for indexBasis = 1:noBasis
                        dRIdx = dRdx(indexBasis);
                        dRIdy = dRdy(indexBasis);
                        
                        id    = (indexBasis-1)*3+1:3*indexBasis;
                        
                        %                         Bmem(:,id)=[dRIdx*a1e1 dRIdx*a1e2 dRIdx*a1e3;
                        %                             dRIa(1,indexBasis)  dRIa(2,indexBasis)  dRIa(3,indexBasis);
                        %                             dRIdy*a2e1 dRIdy*a2e2 dRIdy*a2e3];
                        
                        % ori
                        Bmem(:,id)=[dRIdx*a1e1 dRIdx*a1e2 dRIdx*a1e3;
                            dRIdy*a2e1 dRIdy*a2e2 dRIdy*a2e3;
                            dRIa(1,indexBasis)  dRIa(2,indexBasis)  dRIa(3,indexBasis)];
                        
                        BI1 =   -d2Rdx(indexBasis)*a3 + 1/norma*(dRIdx*cross(a11,a2) + dRIdy*cross(a1,a11) + ...
                            dot(a3,a11)*(dRIdx*cross(a2,a3) + dRIdy*cross(a3,a1)));
                        
                        BI2 = -  d2Rdy(indexBasis)*a3 + 1/norma*(dRIdx*cross(a22,a2) + dRIdy*cross(a1,a22) + ...
                            dot(a3,a22)*(dRIdx*cross(a2,a3) + dRIdy*cross(a3,a1)));
                        
                        BI3 = -  d2Rdxdy(indexBasis)*a3 + 1/norma*(dRIdx*cross(a12,a2) + dRIdy*cross(a1,a12) + ...
                            dot(a3,a12)*(dRIdx*cross(a2,a3) + dRIdy*cross(a3,a1)));
                        
                        %Bben(:,id)=[BI1;2*BI3;BI2];
                        
                        % ori
                        Bben(:,id)=[BI1;BI2;2*BI3];
                    end
%                     disp('size bmem bben')
%                     size(Bmem)
%                     size(Bben)
%                     pause
%                     
                    %                     Bmem
                    %                     pause
                    
                    %                    memStiff * Bmem' * Cmat * Bmem * scalefac * gauss_weight_x(ii).*gauss_weight_y(jj).*J1
%                     Bmem
%                     Bben
%                     pause
                    
                    localstiff = localstiff + ...
                        memStiff * Bmem' * Cmat * Bmem * scalefac * gauss_weight_x(ii).*gauss_weight_y(jj).*J1 + ...
                        benStiff * Bben' * Cmat * Bben * scalefac * gauss_weight_x(ii).*gauss_weight_y(jj).*J1;%
                    
                end
            end
            
            stiff(tsctrx, tsctrx) = stiff(tsctrx, tsctrx) + localstiff;
        end
    end
end
