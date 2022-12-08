function [ basisFun ] = saveBasisShell_plotSol( PHTelem,GIFTmesh, p, q )
%save the basis functions for plot sol
%use with plot solCahnHilliardShell2

%Gauss points
numPts=21;

numPatches = length(PHTelem);
numElements = length(PHTelem{1}); %maybe use the maximum number of elements?
sR = cell(numPatches, numElements, numPts, numPts);
coords = cell(numPatches, numElements, numPts, numPts);

fudge = 0;
uref = linspace(-1+fudge,1-fudge,numPts);
vref = linspace(-1+fudge,1-fudge,numPts);

%1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u,~] = bernstein_basis(uref,p);
[B_v,~] = bernstein_basis(vref,q);
Buv = zeros(numPts, numPts, (p+1)*(q+1));

%the derivatives of the 2D Bernstein polynomials at Gauss points on the
%master element
basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        Buv(:,:,basisCounter) = B_u(:,i)*B_v(:,j)';
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
            for jj=1:numPts
                for ii=1:numPts
                    %evaluate the derivatives of the mapping from parameter
                    %space to physical space
                    
                    sR{patchIndex,elemIndex,ii,jj} = (PHTelem{patchIndex}(elemIndex).modifiedC)*squeeze(Buv(ii,jj,:));      
                    [coord] = paramMapShell(GIFTmesh{patchIndex},uref(ii),vref(jj), xmin, ymin, xmax, ymax);
                    coords{patchIndex,elemIndex,ii,jj} =coord;
                    
                end
            end
        end
    end
end

basisFun.R = sR;
basisFun.coords = coords;



end

