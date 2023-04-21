function plotExactSolBiharmonic1D(PHTelem,GIFTmesh)
u=@(x) -(x*(x-1))^2;
%figure
numPts = 21; %number of plot points to use on each edge
fudge = 0;
uref = linspace(-1+fudge,1-fudge,numPts);

%1D bernstein polynomials evaluated at the Gauss points on the master element
%[B_u,~] = bernstein_basis(uref,p);

for patchIndex = 1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            xmin = PHTelem{patchIndex}(i).vertex(1);
            xmax = PHTelem{patchIndex}(i).vertex(2);

           % globalNodes=PHTelem{patchIndex}(i).nodesGlobal;
            
%            tempDisp=sol(globalNodes);
            coords=zeros(1,numPts);
      %      disp=zeros(1,numPts);
            exactSol=zeros(1,numPts);
            
            for ii=1:numPts
                %evaluate the derivatives of the mapping from parameter
                %space to physical space
                
                [coord,~,~,~] = paramMap1D(GIFTmesh{patchIndex},uref(ii), xmin, xmax );
                coords(ii)=coord;
           %     R = (PHTelem{patchIndex}(i).modifiedC)*(B_u(ii,:)');
                %disp(ii)=R'*tempDisp;
                exactSol(ii)=u(coord);
            end
            
            %   maxZ = max([maxZ; disp(:)]);
            %plot(coords,exactSol,'color','r','LineWidth',2)
            plot(coords,exactSol,'--r','LineWidth',2)
            %surf(coords(:,:,1), coords(:,:,2),disp,'EdgeColor', 'none')
            hold on
            
            
            
        end
    end
end

end

