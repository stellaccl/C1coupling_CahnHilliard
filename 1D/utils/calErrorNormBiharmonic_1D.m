function [l2relerr] =calErrorNormBiharmonic_1D( PHTelem, GIFTmesh, p, sol0)
%cal error biharmonic 1D
u=@(x) -(x*(x-1))^2;  % exact sol

numGaussX = p+1; %number of plot points to use on each edge
[gwx, gpx]=quadrature(numGaussX, 'GAUSS', 1);
%fudge = 0;
%uref = linspace(-1+fudge,1-fudge,numGaussX);
gpx=gpx';
%1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u,~] = bernstein_basis(gpx,p);

l2norm = 0;
l2relerr = 0;

for patchIndex = 1:length(PHTelem)
    
    for i=1:length(PHTelem{patchIndex})
        
        if isempty(PHTelem{patchIndex}(i).children)
            xmin = PHTelem{patchIndex}(i).vertex(1);
            xmax = PHTelem{patchIndex}(i).vertex(2);
            
            scalefac = (xmax - xmin)/2;
            globalNodes=PHTelem{patchIndex}(i).nodesGlobal;
            
            tempDisp=sol0(globalNodes);
            %             coords=zeros(1,numPts);
            %             disp=zeros(1,numPts);
            %loop over the ngauss_x x ngauss_y gauss points on each element
            
            for ii=1:numGaussX
                %evaluate the derivatives of the mapping from parameter
                %space to physical space
                
                [coord,dxdxi,~,~] = paramMap1D(GIFTmesh{patchIndex},gpx(ii), xmin, xmax );

                J = det(dxdxi);
                R = (PHTelem{patchIndex}(i).modifiedC)*(B_u(ii,:)');
                exactSol=u(coord);
                sol=R'*tempDisp;
                
                l2norm = l2norm + (exactSol^2)*gwx(ii)*scalefac*J;  
                l2relerr = l2relerr + ((exactSol-sol)^2)*gwx(ii)*scalefac*J;
                
            end

        end
    end
end

l2relerr = sqrt(l2relerr)/sqrt(l2norm);

end

