function  plotBasis1D_modifiedC_paper(PHUTelem,GIFTmesh,p,indexBasis,c)

numPts = 1001; %number of plot points to use on each edge
fudge = 0;
uref = linspace(-1+fudge,1-fudge,numPts);
colorArray = {'blue', 'red', 'green', 'cyan', 'magenta', 'black','blue', 'red', 'green'};

%1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u,dB_u] = bernstein_basis(uref,p);

R = zeros(numPts, (p+1));
dRdu = zeros(numPts, (p+1));

%the derivatives of the 2D Bernstein polynomials at Gauss points on the
%master element
basisCounter = 0;
for i=1:p+1
    basisCounter = basisCounter + 1;
    R(:,basisCounter) = B_u(:,i);
    dRdu(:,basisCounter) = dB_u(:,i);
end

for indexPatch = 1:length(PHUTelem)
    for indexElem = 1:length(PHUTelem{indexPatch})
        
        if isempty(PHUTelem{indexPatch}(indexElem).children)
            localIndex = (PHUTelem{indexPatch}(indexElem).nodesGlobal==indexBasis);
            
            xmin = PHUTelem{indexPatch}(indexElem).vertex(1);
            xmax = PHUTelem{indexPatch}(indexElem).vertex(2);
            
            if any(localIndex)
                
                cR=zeros(numPts,1);
                dRglob=zeros(numPts,1);
                coords=zeros(numPts,1);
              
                for ii=1:numPts
                    
                    cR(ii) = (PHUTelem{indexPatch}(indexElem).modifiedC(localIndex,:))*R(ii,:)';
                    
                    dR = (PHUTelem{indexPatch}(indexElem).modifiedC(localIndex,:))*dRdu(ii,:)';
                    
                    [coords(ii),dxdxi] = paramMap1D(GIFTmesh{indexPatch}, uref(ii), xmin, xmax );
                    
                    dRglob(ii)=dxdxi\dR;
                    
                    %                     hold on
                    %                     plot(coords,cR,'.','color','r')
                    %                     hold on
                    %                     plot(coords,dRglob,'.','color','b')
                    %                     drawnow
                end
                
                %                     hold on
                %                     plot(coords,cR,'.','color',colorArray{indexBasis})
                %                     hold on
                %                     plot(coords,dRglob,'--','color',colorArray{indexBasis})
                %                     drawnow

                hold on
               % plot(coords,cR,'.','color',c,'MarkerSize',4)
                plot(coords,cR,'color',c,'LineWidth',2)
%                 hold on
%                 plot(coords,dRglob,'--','color',c)
                drawnow
            end
            
        end
    end
end


end

