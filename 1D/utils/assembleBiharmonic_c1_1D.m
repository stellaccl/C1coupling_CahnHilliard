function [stiff, rhs ] = assembleBiharmonic_c1_1D(PHTelem,GIFTmesh, sizeBasis, p )

%Gauss points
ngauss_x = p+1;
[gauss_weight_x, gauss_coord_x] = quadrature( ngauss_x, 'GAUSS', 1 );

%take the transpose so that they are in the format expected by
%bernstein_basis
gauss_coord_x = gauss_coord_x';

%1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u, dB_u,ddB_u] = bernstein_basis(gauss_coord_x,p);

stiff = sparse(sizeBasis,sizeBasis);
rhs= zeros(sizeBasis,1);

for patchIndex = 1:length(PHTelem)
    for elemIndex=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(elemIndex).children)
            xmin = PHTelem{patchIndex}(elemIndex).vertex(1);
            xmax = PHTelem{patchIndex}(elemIndex).vertex(2);
            
            scalefac = (xmax - xmin)/2;
            nument = size(PHTelem{patchIndex}(elemIndex).modifiedC,1);
            scrtx = PHTelem{patchIndex}(elemIndex).nodesGlobal(1:nument);
            
            localStiff = zeros(nument, nument);
            localRhs = zeros(nument, 1);
            
            
            for ii=1:ngauss_x
                
                R = (PHTelem{patchIndex}(elemIndex).modifiedC)*(B_u(ii,:)');
                dRdx = (PHTelem{patchIndex}(elemIndex).modifiedC)*(dB_u(ii,:)');
                d2Rdx = (PHTelem{patchIndex}(elemIndex).modifiedC)*(ddB_u(ii,:)');
                
                dRdx = dRdx*2/(xmax-xmin);
                d2Rdx = d2Rdx*(2/(xmax-xmin))^2;
                
                [coords,dxdxi,d2xdxi2,dxdxi2] = paramMap1D(GIFTmesh{patchIndex}, gauss_coord_x(ii), xmin, xmax );
        
                dR = dxdxi\dRdx';
                ddR=d2Rdx';
                ddR = dxdxi2\(ddR - d2xdxi2*dR);
                J = abs(det(dxdxi));
                
                tempValue= ddR'*ddR;
                localStiff=localStiff+tempValue*scalefac.*gauss_weight_x(ii).*J;
                
                [d4x] = computeDerivatives1DBiharmonic(coords);
                
                f=d4x;
               
                localRhsTemp=R*f;
                
                localRhsTemp=localRhsTemp.*scalefac.*gauss_weight_x(ii).*J;
                
                localRhs = localRhs + localRhsTemp;
                
            end
            
            stiff(scrtx,scrtx) = stiff(scrtx,scrtx) + localStiff;
            rhs(scrtx) =  rhs(scrtx) + localRhs;
        end
    end
end

end