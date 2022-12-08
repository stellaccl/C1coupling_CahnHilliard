function [l2relerr] =calcErrorNormsBiharmonic_squareWith2Holes_degElev(PHTelem, GIFTmesh,  p, q, sol0 )

r_small=sqrt(2)/2;
r_big=1.2;

u=@(x,y)((x-4)^2)*((x+8)^2)*((y-4)^2)*((y+4)^2)*(((x+4)^2+y^2-r_small^2)^2)*((x^2+y^2-r_big^2)^2);


numGaussX = p+2;
numGaussY = q+2;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%
[B_uDE,~,~] = bernstein_basis(gpx,p+1);
[B_vDE,~,~] = bernstein_basis(gpy,q+1);

BuvDE = zeros(numGaussX, numGaussY, (p+2)*(q+2));

basisCounter = 0;
for j=1:q+2
    for i=1:p+2
        basisCounter = basisCounter + 1;
        BuvDE(:,:,basisCounter) = B_uDE(:,i)*B_vDE(:,j)';
      
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
            nument = size(PHTelem{patchIndex}(i).modifiedC,1);
              
            for jj=1:numGaussY
                for ii=1:numGaussX
                    
                    [ coord, dxdxi]  = paramMap( GIFTmesh{patchIndex}, gpx(ii), gpy(jj), xmin, ymin, xmax, ymax);
                   
                     R= zeros(nument,1);
                    
                    for t=1:nument
                        if PHTelem{patchIndex}(i).polyDegree(t)==p 
                            R(t) = (PHTelem{patchIndex}(i).modifiedC(t,1:(p+1)*(q+1)))*squeeze(Buv(ii,jj,:));      
                        else  
                            R(t) = (PHTelem{patchIndex}(i).modifiedC(t,:))*squeeze(BuvDE(ii,jj,:));       
                        end
                    end

                    
                    J = det(dxdxi);
                   
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

