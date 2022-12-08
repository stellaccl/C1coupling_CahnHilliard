function PlotDispShellC1( PHTelem, GIFTmesh,  p, q, sol)
%plots the X,Y,Z displacements for the thin shell problems
%uses GIFT mapping
%supports multipatches
figure
numPts = 2; %number of plot points to use on each edge
fudge = 0;
uref = linspace(-1+fudge,1-fudge,numPts);
vref = linspace(-1+fudge,1-fudge,numPts);
minZ = Inf;

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

%factor = 1e6;
factor = 1e3;
indexMat = [1,2;4,3]; %index of corners coresponding to (ii,jj)

for patchIndex = 1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            xmin = PHTelem{patchIndex}(i).vertex(1);
            xmax = PHTelem{patchIndex}(i).vertex(3);
            ymin = PHTelem{patchIndex}(i).vertex(2);
            ymax = PHTelem{patchIndex}(i).vertex(4);
            
            sctr=PHTelem{patchIndex}(i).nodesGlobal;
            sctr_x = 3*sctr-2;
            sctr_y = 3*sctr-1;
            sctr_z = 3*sctr;
          %  tsctrx = reshape([sctr_x; sctr_y; sctr_z],1,3*nument);
            
            tempDisp_x=sol(sctr_x);
            tempDisp_y=sol(sctr_y);
            tempDisp_z=sol(sctr_z);
            
            coords_x=zeros(numPts*numPts,1);
            coords_y=zeros(numPts*numPts,1);
            coords_z=zeros(numPts*numPts,1);
            
            disp_x=zeros(numPts*numPts,1);
            disp_y=zeros(numPts*numPts,1);
            disp_z=zeros(numPts*numPts,1);
            %loop over the corners of each element
            
            for jj=1:numPts
                for ii=1:numPts
                    %evaluate the derivatives of the mapping from parameter
                    %space to physical space
                    index = indexMat(ii,jj);
               
                    [ coord] = paramMapShell( GIFTmesh{patchIndex}, uref(ii), vref(jj), xmin, ymin, xmax, ymax);
                    coords_x(index) = coord(1);
                    coords_y(index) = coord(2);
                    coords_z(index) = coord(3);    
                    R = (PHTelem{patchIndex}(i).modifiedC)*squeeze(Buv(ii,jj,:));
                    disp_x(index)=R'*tempDisp_x;
                    disp_y(index)=R'*tempDisp_y;
                    disp_z(index)=R'*tempDisp_z;
                    
                end
            end
            %minZ = min([minZ; disp(:)]);
            fill3(coords_x+factor*disp_x, coords_y+factor*disp_y,coords_z+factor*disp_z, disp_y,'EdgeColor','none')
            hold on
           
        end
    end
end
