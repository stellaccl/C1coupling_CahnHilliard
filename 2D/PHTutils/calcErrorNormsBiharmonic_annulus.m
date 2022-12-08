function [l2relerr] =calcErrorNormsBiharmonic_annulus(PHTelem, GIFTmesh,  p, q, sol0 )
r_in=0.2;
r_out=1;
u=@(x,y) ((x^2+y^2-r_in^2)^2)*((x^2+y^2-r_out^2)^2);


disp('enter cal error')
numGaussX = p+1;
numGaussY = q+1;

[gwx, gpx]=quadrature(numGaussX, 'GAUSS', 1);
[gwy, gpy]=quadrature(numGaussY, 'GAUSS', 1);

gpx=gpx';
gpy=gpy';

l2norm = 0;
l2relerr = 0;

%define the 2D Bernstein polynomials
[B_u, ~] = bernstein_basis(gpx,p);
[B_v, ~] = bernstein_basis(gpy,q);
Buv = zeros(numGaussX, numGaussY,(p+1)*(q+1));

basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        Buv(:,:,basisCounter) = B_u(:,i)*B_v(:,j)';
    end
end

for patchIndex = 1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            xmin = PHTelem{patchIndex}(i).vertex(1);
            xmax = PHTelem{patchIndex}(i).vertex(3);
            ymin = PHTelem{patchIndex}(i).vertex(2);
            ymax = PHTelem{patchIndex}(i).vertex(4);
            
            scalefac = (xmax - xmin)*(ymax - ymin)/4;
            globalNodes=PHTelem{patchIndex}(i).nodesGlobal;
            tempSol=sol0(globalNodes);
            
            for jj=1:numGaussY
                for ii=1:numGaussX
                    
                    [ coord, dxdxi]  = paramMap( GIFTmesh{patchIndex}, gpx(ii), gpy(jj), xmin, ymin, xmax, ymax);
                    
                    J = det(dxdxi);
                    R = (PHTelem{patchIndex}(i).modifiedC)*squeeze(Buv(ii,jj,:));
                    sol=R'*tempSol;
                    %exactSol = coord(1)^2*(1-coord(1))^2*coord(2)^2*(1-coord(2))^2;
                    exactSol =u(coord(1),coord(2));
                    l2norm = l2norm + (exactSol^2)*gwx(ii)*gwy(jj)*scalefac*J;
                    l2relerr = l2relerr + ((exactSol-sol)^2)*gwx(ii)*gwy(jj)*scalefac*J;
                    
                end
            end
            
        end
    end
end
%relative error
l2relerr = sqrt(l2relerr)/sqrt(l2norm);

end

