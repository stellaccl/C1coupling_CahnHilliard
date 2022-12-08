function  PlotBasisShell_multiElemGIFTmesh( PHUTelem,GIFTmesh,p,q,type2Basis)

figure
numPts = 11; %number of plot points to use on each edge
fudge = 0;
uref = linspace(-1+fudge,1-fudge,numPts);
vref = linspace(-1+fudge,1-fudge,numPts);

%1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u,dB_u,ddB_u] = bernstein_basis(uref,p);
[B_v,dB_v,ddB_v] = bernstein_basis(vref,q);

Buv = zeros(numPts, numPts, (p+1)*(q+1));
dBdu = zeros(numPts, numPts, (p+1)*(q+1));
dBdv = zeros(numPts, numPts, (p+1)*(q+1));
d2Bdu = zeros(numPts, numPts, (p+1)*(q+1));
d2Bdv = zeros(numPts, numPts, (p+1)*(q+1));
d2Bdudv = zeros(numPts, numPts, (p+1)*(q+1));
%the derivatives of the 2D Bernstein polynomials at Gauss points on the
%master element
basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        Buv(:,:,basisCounter) = B_u(:,i)*B_v(:,j)';
        dBdu(:,:,basisCounter) = dB_u(:,i)*B_v(:,j)';
        dBdv(:,:,basisCounter) = B_u(:,i)*dB_v(:,j)';
        d2Bdu(:,:,basisCounter) = ddB_u(:,i)*B_v(:,j)';
        d2Bdv(:,:,basisCounter) = B_u(:,i)*ddB_v(:,j)';
        d2Bdudv(:,:,basisCounter) = dB_u(:,i)*dB_v(:,j)';
    end
end

R = zeros(numPts, numPts);
derivativesX = zeros(numPts, numPts);
derivativesY = zeros(numPts, numPts);
derivativesZ = zeros(numPts, numPts);

coords=zeros(numPts,numPts,3);


numPatches = length(PHUTelem);
for indexBasis=1:length(type2Basis)
    for indexPatch = 1:numPatches
        for indexElem = 1:length(PHUTelem{indexPatch})
            
            
            up=PHUTelem{indexPatch}(indexElem).neighbor_up;
            down=PHUTelem{indexPatch}(indexElem).neighbor_down;
            left=PHUTelem{indexPatch}(indexElem).neighbor_left;
            right=PHUTelem{indexPatch}(indexElem).neighbor_right;
            
            if isempty(PHUTelem{indexPatch}(indexElem).children)
                if ismember(type2Basis(indexBasis),PHUTelem{indexPatch}(indexElem).nodesGlobal)
                    
                    localIndex = (PHUTelem{indexPatch}(indexElem).nodesGlobal==type2Basis(indexBasis));
                    
                    xmin = PHUTelem{indexPatch}(indexElem).vertex(1);
                    xmax = PHUTelem{indexPatch}(indexElem).vertex(3);
                    ymin = PHUTelem{indexPatch}(indexElem).vertex(2);
                    ymax = PHUTelem{indexPatch}(indexElem).vertex(4);
                    %
                    
                    for jj=1:numPts
                        for ii=1:numPts
                            
                            % [coord, dxdxi, d2xdxi2] =paramMapSemiSphere3(PHUTelem, indexPatch,indexElem,uref(ii),vref(jj));
                            % [coord, dxdxi, d2xdxi2] =paramMapSphere4(PHUTelem, indexPatch,indexElem,uref(ii),vref(jj));
                            [ coord, dxdxi] = paramMapShell_multiElemGIFTmesh(indexElem,GIFTmesh{indexPatch}, uref(ii),vref(jj), xmin, ymin, xmax, ymax);
                            coords(jj,ii,:)=[coord(1),coord(2),coord(3)];
                            R(jj,ii) = (PHUTelem{indexPatch}(indexElem).modifiedC(localIndex,:))*squeeze( Buv(ii,jj,:));
                            
                            dRdx = (PHUTelem{indexPatch}(indexElem).modifiedC(localIndex,:))*squeeze(dBdu(ii,jj,:));
                            dRdy = (PHUTelem{indexPatch}(indexElem).modifiedC(localIndex,:))*squeeze(dBdv(ii,jj,:));
                            
                            dRdx = dRdx*2/(xmax-xmin);
                            dRdy = dRdy*2/(ymax-ymin);
                            
                            J=dxdxi';
                            G=J'*J;
                            dR=J*(inv(G))'*([dRdx,dRdy]');
                            %
                            derivativesX(jj,ii)=dR(1);
                            derivativesY(jj,ii)=dR(2);
                            derivativesZ(jj,ii)=dR(3);
                            
                        end
                    end
                    
                                       surf(coords(:,:,1), coords(:,:,2),coords(:,:,3),R,'FaceColor','interp','EdgeColor','none')
                                        alpha(0.75)
                                        title(['basis function',num2str(type2Basis(indexBasis))])
                                        hold on
                                        axis equal
                    
                    
                    %                     subplot(2,2,2)
%                     surf(coords(:,:,1), coords(:,:,2),coords(:,:,3), derivativesX,'FaceColor','interp','EdgeColor','none')
%                     title('derivatives x')
%                     hold on
%                     axis equal
                    %
                    %
                    %                     subplot(2,2,3)
%                                         surf(coords(:,:,1), coords(:,:,2),coords(:,:,3), derivativesY,'FaceColor','interp','EdgeColor','none')
%                                         title('derivatives y')
%                                         hold on
%                                         axis equal
                    % %
                    %
                    %                     subplot(2,2,4)
%                                         surf(coords(:,:,1), coords(:,:,2),coords(:,:,3), derivativesZ,'FaceColor','interp','EdgeColor','none')
%                                         title('derivatives z')
%                                         hold on
%                                         axis equal
                    
                    
                    %                     plot patch boundary
                    coord_south = zeros(numPts,3);
                    coord_east = zeros(numPts,3);
                    coord_north = zeros(numPts,3);
                    coord_west = zeros(numPts,3);
                    
                    for j=1:numPts
                        
                        coord_south(j,:) = paramMapShell_multiElemGIFTmesh(indexElem,GIFTmesh{indexPatch}, uref(j), -1,  xmin, ymin, xmax, ymax );
                        coord_east(j,:) = paramMapShell_multiElemGIFTmesh(indexElem,GIFTmesh{indexPatch}, 1, vref(j), xmin, ymin, xmax, ymax);
                        coord_north(j,:) = paramMapShell_multiElemGIFTmesh(indexElem,GIFTmesh{indexPatch}, -uref(j), 1,  xmin, ymin, xmax, ymax );
                        coord_west(j,:) = paramMapShell_multiElemGIFTmesh(indexElem,GIFTmesh{indexPatch}, -1, -vref(j),  xmin, ymin, xmax, ymax );
                    end
                    
                    if isempty(up)
                        
                        line( coord_north (:,1), coord_north(:,2), coord_north(:,3), 'Color', 'k')
                    end
                    
                    if isempty(right)
                        
                        line( coord_east (:,1), coord_east(:,2), coord_east(:,3), 'Color', 'k')
                    end
                    
                    
                    if isempty(down)
                        
                        line( coord_south (:,1), coord_south(:,2), coord_south(:,3), 'Color', 'k')
                    end
                    
                    
                    if isempty(left)
                        
                        line( coord_west(:,1), coord_west(:,2), coord_west(:,3), 'Color', 'k')
                    end
                    %subplot(2,2,1)
                    
                    
                    
                    
                end
            end
        end
    end
    
    %plotPHTMeshSphere_new2( PHUTelem)
   % close all
end










end
