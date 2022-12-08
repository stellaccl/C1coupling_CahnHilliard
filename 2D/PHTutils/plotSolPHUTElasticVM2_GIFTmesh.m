function  plotSolPHUTElasticVM2_GIFTmesh( sol0, GIFTmesh,PHUTelem, p, q, Cmat)
%plot subroutine for spannerElasticity.m
%plots the von Misses stresses
%for use with patchTestLshapeQuartic...

numPts =2;

%plots the deformed shape + stresses
fudge=0;
xi = linspace(-1+fudge,1-fudge,numPts);
eta = linspace(-1+fudge,1-fudge,numPts);

%calculate the number of actual elements (i.e., non-refined, without children)
numElem = 0;
for indexPatch = 1:length(PHUTelem)
    for i=1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(i).children)
            numElem = numElem+1;
        end
    end
end

%Displaying the displacements
element4 = zeros(numElem, 4);
physcoord = zeros(4*numElem, 2);
dispcoord = zeros(4*numElem, 2);
straincoord = zeros(4*numElem, 3);
sigmacoord = zeros(4*numElem, 3);

%define the 2D Bernstein polynomials
[B_u, dB_u] = bernstein_basis(xi,p);
[B_v, dB_v] = bernstein_basis(eta,q);

Buv = zeros(numPts, numPts, numPts*numPts);
dBdu = zeros(numPts, numPts, numPts*numPts);
dBdv = zeros(numPts, numPts, numPts*numPts);

basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        Buv(:,:,basisCounter) = B_u(:,i)*B_v(:,j)';
        dBdu(:,:,basisCounter) = dB_u(:,i)*B_v(:,j)';
        dBdv(:,:,basisCounter) = B_u(:,i)*dB_v(:,j)';
    end
end

elementCounter = 0;
for indexPatch = 1:length(PHUTelem)
    for i=1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(i).children) 

            elementCounter = elementCounter + 1;
            element4(elementCounter, :) = (elementCounter-1)*4+1:(elementCounter-1)*4+4;
            
            coord = cell(numPts,numPts);
            
            nument = size(PHUTelem{indexPatch}(i).modifiedC,1); %number of basis functions with support on current knotspan
            
            %initialize matrices to store the displacement, strain and stress
            %values at each plot point
            dispmatx = zeros(numPts,numPts);
            dispmaty = zeros(numPts,numPts);
            strain11 = zeros(numPts,numPts);
            strain12 = zeros(numPts,numPts);
            strain22 = zeros(numPts,numPts);
            stressvect = cell(numPts,numPts);

            scrt = PHUTelem{indexPatch}(i).nodesGlobal(1:nument);
            scrt_x = 2*scrt-1;
            scrt_y = 2*scrt;
            dscrtx = reshape([2*scrt-1; 2*scrt],1,2*nument);
            
            xmin = PHUTelem{indexPatch}(i).vertex(1);
            xmax = PHUTelem{indexPatch}(i).vertex(3);
            ymin = PHUTelem{indexPatch}(i).vertex(2);
            ymax = PHUTelem{indexPatch}(i).vertex(4);

            for jj=1:numPts
                for ii=1:numPts
                    
                    %compute the mapping from reference space to physical space
                    
                    %coord_x = xi(ii)*(xmax-xmin)/2+(xmax+xmin)/2;
                    %coord_y = eta(jj)*(ymax-ymin)/2+(ymax+ymin)/2;
                    %coord_pt = [coord_x, coord_y];
                    
                    %evaluate the basis functions
                    
                    cR = PHUTelem{indexPatch}(i).modifiedC * squeeze(Buv(ii,jj,:));
                    cdRdx = PHUTelem{indexPatch}(i).modifiedC * squeeze(dBdu(ii,jj,:))*2/(xmax-xmin);
                    cdRdy = PHUTelem{indexPatch}(i).modifiedC * squeeze(dBdv(ii,jj,:))*2/(ymax-ymin);
                    
                    u_hat = xi(ii);
                    v_hat = eta(jj);
                                       
                   %[coord_pt,~,~, dxdxi] = paramMapPHUT(u_hat, v_hat, verts);
                   
                   [coord_pt, dxdxi] = paramMap( GIFTmesh{indexPatch}, u_hat, v_hat, xmin, ymin, xmax, ymax);
              
                   
                    %dxdxi = dxdxi';
                    coord{jj,ii} = [coord_pt(1), coord_pt(2)];
                    
                    dRloc = [cdRdx';cdRdy'];
                    
                    if abs(det(dxdxi))>1e-4
                        dR =  dxdxi\dRloc;
                    else
                        dR = pinv(dxdxi)*dRloc;
                    end
                    
                    B = zeros(2*nument,3);
                    B(1:2:2*nument-1,1) = dR(1,:);
                    B(2:2:2*nument,2) = dR(2,:);
                    B(1:2:2*nument-1,3) = dR(2,:);
                    B(2:2:2*nument,3) = dR(1,:);
                    
                    %calculate displacement values
                    dispmatx(jj,ii) = dispmatx(jj,ii) + cR'*sol0(scrt_x);
                    dispmaty(jj,ii) = dispmaty(jj,ii) + cR'*sol0(scrt_y);
                    
                    %calculate the stress values
                    stressvect{jj,ii} = Cmat*B'*sol0(dscrtx);
                    
                    
                end
            end
            
            physcoord((elementCounter-1)*4+1, :) = coord{1,1};
            physcoord((elementCounter-1)*4+2, :) = coord{1,end};
            physcoord((elementCounter-1)*4+3, :) = coord{end,end};
            physcoord((elementCounter-1)*4+4, :) = coord{end,1};
            
            dispcoord((elementCounter-1)*4+1, :) = [dispmatx(1,1) dispmaty(1,1)];
            dispcoord((elementCounter-1)*4+2, :) = [dispmatx(1,2) dispmaty(1,2)];
            dispcoord((elementCounter-1)*4+3, :) = [dispmatx(2,2) dispmaty(2,2)];
            dispcoord((elementCounter-1)*4+4, :) = [dispmatx(2,1) dispmaty(2,1)];
            straincoord((elementCounter-1)*4+1, :) = [strain11(1,1) strain12(1,1) strain22(1,1)];
            straincoord((elementCounter-1)*4+2, :) = [strain11(1,2) strain12(1,2) strain22(1,2)];
            straincoord((elementCounter-1)*4+3, :) = [strain11(2,2) strain12(2,2) strain22(2,2)];
            straincoord((elementCounter-1)*4+4, :) = [strain11(2,1) strain12(2,1) strain22(2,1)];
            
            
            sigmacoord((elementCounter-1)*4+1, :) = stressvect{1,1}';
            sigmacoord((elementCounter-1)*4+2, :) = stressvect{1,2}';
            sigmacoord((elementCounter-1)*4+3, :) = stressvect{2,2}';
            sigmacoord((elementCounter-1)*4+4, :) = stressvect{2,1}';
            
            hold on
        end
    end
end
disp('Trisurfing...')
factor = 1;

figure
%trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sqrt(sigmacoord(:,1).^2+sigmacoord(:,2).^2-sigmacoord(:,1).*sigmacoord(:,2)+3*sigmacoord(:,3).^2), 'EdgeColor','none','facecolor','interp')
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sqrt(sigmacoord(:,1).^2+sigmacoord(:,2).^2-sigmacoord(:,1).*sigmacoord(:,2)), 'EdgeColor','none','facecolor','interp')
view(0,90)
title('Displacements and \sigma_{VM}')
colorbar('vert')
drawnow


figure
factor = 0;
%trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sqrt(sigmacoord(:,1).^2+sigmacoord(:,2).^2-sigmacoord(:,1).*sigmacoord(:,2)+3*sigmacoord(:,3).^2), 'EdgeColor','none','facecolor','interp')
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), dispcoord(:,1), 'EdgeColor','none','facecolor','interp')
view(0,90)
title('Displacements in u direction')
colorbar('vert')
drawnow


figure
factor = 0;
%trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sqrt(sigmacoord(:,1).^2+sigmacoord(:,2).^2-sigmacoord(:,1).*sigmacoord(:,2)+3*sigmacoord(:,3).^2), 'EdgeColor','none','facecolor','interp')
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sqrt(dispcoord(:,1).^2+dispcoord(:,2).^2), 'EdgeColor','none','facecolor','interp')
view(0,90)
title('Displacement magnitude')
colorbar('vert')
drawnow

