function [] = PlotExactSol_annulusBiharmonic(PHTelem, GIFTmesh,  p, q)
r_in=0.2;
r_out=1;
u=@(x,y) ((x^2+y^2-r_in^2)^2)*((x^2+y^2-r_out^2)^2); 
%u=((x^2+y^2-r_in^2)^2)*((x^2+y^2-r_out^2)^2)

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
            
            globalNodes=PHTelem{patchIndex}(i).nodesGlobal;
            %tempDisp=sol(globalNodes);
            coords=zeros(numPts,numPts,2);
            exactSol=zeros(numPts,numPts);
            %loop over the ngauss_x x ngauss_y gauss points on each element
            for jj=1:numPts
                for ii=1:numPts
                    %evaluate the derivatives of the mapping from parameter
                    %space to physical space
                    
                    [ coord, ~] = paramMapPlate( GIFTmesh{patchIndex},uref(ii), vref(jj), xmin, ymin, xmax, ymax);
                    %J = det(dxdxi);
                    coords(ii,jj,:)=[coord(1),coord(2)];
                    %R = (PHTelem{patchIndex}(i).modifiedC)*squeeze(Buv(ii,jj,:));
                    %disp=R'*tempDisp;
                    
                    tempExactSol = u(coord(1),coord(2));
                    %tempExactSol = coord(1)^2*(1-coord(1))^2*coord(2)^2*(1-coord(2))^2;
                    %diff(ii,jj)=abs(exactSol-disp);
                    exactSol(ii,jj)=tempExactSol;
                end
            end
            maxZ = max([maxZ; exactSol(:)]);
            surf(coords(:,:,1), coords(:,:,2),exactSol,'EdgeColor', 'none')
            hold on
            
        end
    end
end

end

