function PlotDisp1D( PHTelem, GIFTmesh,  p, sol)
%plots the Z displacement for the plate problem and calculates the minimum
%displacement from the plotting data
%uses GIFT mapping
%supports multipatches
%use modified C

figure
numPts = 21; %number of plot points to use on each edge
fudge = 0;
uref = linspace(-1+fudge,1-fudge,numPts);

%1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u,~] = bernstein_basis(uref,p);

for patchIndex = 1:length(PHTelem)
    
    for i=1:length(PHTelem{patchIndex})
        
        if isempty(PHTelem{patchIndex}(i).children)
            xmin = PHTelem{patchIndex}(i).vertex(1);
            xmax = PHTelem{patchIndex}(i).vertex(2);
            
            
            globalNodes=PHTelem{patchIndex}(i).nodesGlobal;
    
            tempDisp=sol(globalNodes);
            coords=zeros(1,numPts);
            disp=zeros(1,numPts);
            %loop over the ngauss_x x ngauss_y gauss points on each element
            
            for ii=1:numPts
                %evaluate the derivatives of the mapping from parameter
                %space to physical space
                
                [coord,~,~,~] = paramMap1D(GIFTmesh{patchIndex},uref(ii), xmin, xmax );
                coords(ii)=coord;
                R = (PHTelem{patchIndex}(i).modifiedC)*(B_u(ii,:)');
                
                disp(ii)=R'*tempDisp;
                
            end
            
            %   maxZ = max([maxZ; disp(:)]);
            plot(coords,disp,'color','k','LineWidth',2)
            %surf(coords(:,:,1), coords(:,:,2),disp,'EdgeColor', 'none')
            hold on
            
        end
    end
end
