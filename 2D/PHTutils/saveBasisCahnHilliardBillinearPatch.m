function [ basisFun ] = saveBasisCahnHilliardBillinearPatch( PHTelem,GIFTmesh, p, q )
%save the basis functions for Cahn Hilliard equation
disp('in saveBasis')
%Gauss points
ngauss_x = p+1;
ngauss_y = q+1;

numPatches = length(PHTelem);
numElements = length(PHTelem{1}); %maybe use the maximum number of elements?
sR = cell(numPatches, numElements, ngauss_x, ngauss_y);
sdR = cell(numPatches, numElements, ngauss_x, ngauss_y);
sdRdx = cell(numPatches, numElements, ngauss_x, ngauss_y);
sdRdy = cell(numPatches, numElements, ngauss_x, ngauss_y);
sd2Rdx = cell(numPatches, numElements, ngauss_x, ngauss_y);                
sd2Rdy = cell(numPatches, numElements, ngauss_x, ngauss_y);
sd2Rdxdy = cell(numPatches, numElements, ngauss_x, ngauss_y);


slaplacian= cell(numPatches, numElements, ngauss_x, ngauss_y);
sdxdxi= cell(numPatches, numElements, ngauss_x, ngauss_y);
sd2xdxi2= cell(numPatches, numElements, ngauss_x, ngauss_y);
sdxdxi2= cell(numPatches, numElements, ngauss_x, ngauss_y);
sJ1= cell(numPatches, numElements, ngauss_x, ngauss_y);

[~, gauss_coord_x] = quadrature( ngauss_x, 'GAUSS', 1 );
[~, gauss_coord_y] = quadrature( ngauss_y, 'GAUSS', 1 );

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

for patchIndex = 1:length(PHTelem)
    for elemIndex=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(elemIndex).children)
            xmin = PHTelem{patchIndex}(elemIndex).vertex(1);
            xmax = PHTelem{patchIndex}(elemIndex).vertex(3);
            ymin = PHTelem{patchIndex}(elemIndex).vertex(2);
            ymax = PHTelem{patchIndex}(elemIndex).vertex(4);
            
            %the jacobian of the transformation from [-1,1]x[-1,1] to
            %[xmin, xmax]x [ymin, ymax]
            
            %elemIndex
            %loop over the ngauss_x x ngauss_y gauss points on each element
            for jj=1:ngauss_y
                for ii=1:ngauss_x
                    %evaluate the derivatives of the mapping from parameter
                    %space to physical space
                    
                    sR{patchIndex,elemIndex,ii,jj} = (PHTelem{patchIndex}(elemIndex).modifiedC)*squeeze(Buv(ii,jj,:));
                    dRdx = (PHTelem{patchIndex}(elemIndex).modifiedC)*squeeze(dBdu(ii,jj,:));
                    dRdy = (PHTelem{patchIndex}(elemIndex).modifiedC)*squeeze(dBdv(ii,jj,:));
                    d2Rdx = (PHTelem{patchIndex}(elemIndex).modifiedC)*squeeze(d2Bdu(ii,jj,:));
                    d2Rdy = (PHTelem{patchIndex}(elemIndex).modifiedC)*squeeze(d2Bdv(ii,jj,:));
                    d2Rdxdy = (PHTelem{patchIndex}(elemIndex).modifiedC)*squeeze(d2Bdudv(ii,jj,:));
                    
                    %[~, dxdxi, d2xdxi2] =  paramMapSemiSphere3(PHTelem, patchIndex,elemIndex, uref, vref);
                    dRdx = dRdx*2/(xmax-xmin);
                    dRdy = dRdy*2/(ymax-ymin);
                    d2Rdx = d2Rdx*(2/(xmax-xmin))^2;
                    d2Rdy = d2Rdy*(2/(ymax-ymin))^2;
                    d2Rdxdy = d2Rdxdy*(2/(xmax-xmin))*(2/(ymax-ymin));

                    [~, dxdxi, d2xdxi2, dxdxi2] = paramMapPlate( GIFTmesh{patchIndex},gauss_coord_x(ii),gauss_coord_y(jj), xmin, ymin, xmax, ymax);

                    J = det(dxdxi);
                    dR = dxdxi\[dRdx'; dRdy'];
                    ddR=[d2Rdx';d2Rdxdy';d2Rdy'];
                    ddR = dxdxi2\(ddR - d2xdxi2*dR);

                    sdR{patchIndex,elemIndex,ii,jj} =dR;
                    sdRdx{patchIndex,elemIndex,ii,jj} =dRdx;
                    sdRdy{patchIndex,elemIndex,ii,jj} =dRdy;
                    
                    sd2Rdx{patchIndex,elemIndex,ii,jj} =ddR(1,:);
                    sd2Rdy{patchIndex,elemIndex,ii,jj} =ddR(3,:);
                    sd2Rdxdy{patchIndex,elemIndex,ii,jj} =ddR(2,:);
                    
                    sdxdxi{patchIndex,elemIndex,ii,jj} = dxdxi;
                    sd2xdxi2{patchIndex,elemIndex,ii,jj} = d2xdxi2;
                    sdxdxi2{patchIndex,elemIndex,ii,jj} =  dxdxi2;
%                    % J=dxdxi';
                   % G=J'*J;
                   % sdR{patchIndex,elemIndex,ii,jj} =J*(inv(G))'*([dRdx,dRdy]');
                  
                    % [laplacian] = computeLaplaceBeltrami(G,dxdxi,d2xdxi2,dRdx,dRdy,d2Rdx,d2Rdy,d2Rdxdy);
                    
                   laplacian=(ddR(1,:)'+ddR(3,:)');
                    
                   slaplacian{patchIndex,elemIndex,ii,jj} = laplacian;
                    
%                     a1    = dxdxi(1,:);
%                     a2    = dxdxi(2,:);
%                     a3    = cross(a1,a2);
%                     norma = norm(a3);
%                    sJ1{patchIndex,elemIndex,ii,jj}    = norma;
                    sJ1{patchIndex,elemIndex,ii,jj}    =  det(J);
                    
                end
            end
        end
    end
end

basisFun.R = sR;
basisFun.dRdx = sdRdx;
basisFun.dRdy = sdRdy;
basisFun.d2Rdx = sd2Rdx;
basisFun.d2Rdy = sd2Rdy;
basisFun.d2Rdxdy = sd2Rdxdy;
basisFun.dxdxi = sdxdxi;
basisFun.d2xdxi2 = sd2xdxi2;
basisFun.dxdxi2 = sdxdxi2;
basisFun.dR = sdR;
basisFun.laplacian = slaplacian;
basisFun.J1 = sJ1;

% save('sdRsaveBasis.mat','sdR')
% disp('finishSaveBasis')
% pause
end

