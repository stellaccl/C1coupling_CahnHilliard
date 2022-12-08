function  plotExactSol_circle_biharmonic( PHTelem, GIFTmesh,  p, q)
%plots the Z displacement for the plate problem and calculates the minimum
%displacement from the plotting data
%uses GIFT mapping
%supports multipatches
%use modified C
r=0.5;
u=@(x,y) (x^2+y^2-r^2)^2;


figure
numPts = 21; %number of plot points to use on each edge
fudge = 0;
uref = linspace(-1+fudge,1-fudge,numPts);
vref = linspace(-1+fudge,1-fudge,numPts);
maxZ = -Inf;

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
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            xmin = PHTelem{patchIndex}(i).vertex(1);
            xmax = PHTelem{patchIndex}(i).vertex(3);
            ymin = PHTelem{patchIndex}(i).vertex(2);
            ymax = PHTelem{patchIndex}(i).vertex(4);
            
            %globalNodes=PHTelem{patchIndex}(i).nodesGlobal;
           % tempDisp=sol(globalNodes);
            coords=zeros(numPts,numPts,2);
            disp=zeros(numPts,numPts);
            %loop over the ngauss_x x ngauss_y gauss points on each element
            for jj=1:numPts
                for ii=1:numPts
                    %evaluate the derivatives of the mapping from parameter
                    %space to physical space
                    
                    [ coord, dxdxi] = paramMapPlate( GIFTmesh{patchIndex},uref(ii), vref(jj), xmin, ymin, xmax, ymax);
                   %  J = abs(det(dxdxi));
                    coords(ii,jj,:)=[coord(1),coord(2)];
                    exactSol =u(coord(1),coord(2));
                    %R = (PHTelem{patchIndex}(i).modifiedC)*squeeze(Buv(ii,jj,:));
                    disp(ii,jj)=exactSol;
                end
            end
          %  maxZ = max([maxZ; disp(:)]);
          
            surf(coords(:,:,1), coords(:,:,2),disp,'EdgeColor', 'none')
            hold on
            
        end
    end
end
