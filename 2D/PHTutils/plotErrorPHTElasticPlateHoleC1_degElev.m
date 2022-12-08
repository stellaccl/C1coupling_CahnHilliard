function  plotErrorPHTElasticPlateHoleC1_degElev( sol0, PHTelem, GIFTmesh, p, q, Cmat, rad, tx )
%plots the error in stresses for the plate with hole example
%supports multipatches

numPatches = length(PHTelem);
numPts = 2;
fudge = 1e-3;
%plots the deformed shape + stresses
xi = linspace(-1+fudge,1-fudge,numPts);
eta = linspace(-1+fudge,1-fudge,numPts);


%calculate the number of actual elements (i.e., non-refined, without children)
numElem = 0;
for patchIndex=1:numPatches
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            numElem = numElem+1;
        end
    end
end


%Displaying the displacements
element4 = zeros(numElem, 4);
physcoord = zeros(4*numElem, 2);
dispcoord = zeros(4*numElem, 2);
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

%define the degree elevated 2D Bernstein polynomials
[B_uDE, dB_uDE] = bernstein_basis(xi,p+1);
[B_vDE, dB_vDE] = bernstein_basis(eta,q+1);

BuvDE = zeros(numPts, numPts, numPts*numPts);
dBduDE = zeros(numPts, numPts, numPts*numPts);
dBdvDE = zeros(numPts, numPts, numPts*numPts);

basisCounter = 0;
for j=1:q+2
    for i=1:p+2
        basisCounter = basisCounter + 1;
        BuvDE(:,:,basisCounter) = B_uDE(:,i)*B_vDE(:,j)';
        dBduDE(:,:,basisCounter) = dB_uDE(:,i)*B_vDE(:,j)';
        dBdvDE(:,:,basisCounter) = B_uDE(:,i)*dB_vDE(:,j)';
    end
end


elementCounter = 0;
for patchIndex=1:numPatches
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            
            elementCounter = elementCounter + 1;
            element4(elementCounter, :) = (elementCounter-1)*4+1:(elementCounter-1)*4+4;
            xmin = PHTelem{patchIndex}(i).vertex(1);
            xmax = PHTelem{patchIndex}(i).vertex(3);
            ymin = PHTelem{patchIndex}(i).vertex(2);
            ymax = PHTelem{patchIndex}(i).vertex(4);
            coord = cell(numPts,numPts);
            
            nument = size(PHTelem{patchIndex}(i).modifiedC,1); %number of basis functions with support on current knotspan
            
            %initialize matrices to store the displacement, strain and stress
            %values at each plot point
            dispmatx = zeros(numPts,numPts);
            dispmaty = zeros(numPts,numPts);
            stressvect = cell(numPts,numPts);
            
            scrt = PHTelem{patchIndex}(i).nodesGlobal(1:nument);
            scrt_x = 2*scrt-1;
            scrt_y = 2*scrt;
            dscrtx = reshape([2*scrt-1; 2*scrt],1,2*nument);
            
            for jj=1:numPts
                for ii=1:numPts
                    
                    %compute the mapping from parameter space to physical space
                    [ coord_pt, dxdxi]  = paramMap( GIFTmesh{patchIndex}, xi(ii), eta(jj), xmin, ymin, xmax, ymax);
                    
                    %evaluate the basis functions
                    cR = zeros(nument,1);
                    cdRdx = zeros(nument,1);
                    cdRdy = zeros(nument,1);
                    for t=1:nument
                        if PHTelem{patchIndex}(i).polyDegree(t)==p
                            cR(t) = PHTelem{patchIndex}(i).modifiedC(t,1:(p+1)*(q+1)) * squeeze(Buv(ii,jj,:));
                            cdRdx(t) = PHTelem{patchIndex}(i).modifiedC(t,1:(p+1)*(q+1)) * squeeze(dBdu(ii,jj,:))*2/(xmax-xmin);
                            cdRdy(t) = PHTelem{patchIndex}(i).modifiedC(t,1:(p+1)*(q+1)) * squeeze(dBdv(ii,jj,:))*2/(ymax-ymin);
                        else
                            cR(t) = PHTelem{patchIndex}(i).modifiedC(t,:) * squeeze(BuvDE(ii,jj,:));
                            cdRdx(t) = PHTelem{patchIndex}(i).modifiedC(t,:) * squeeze(dBduDE(ii,jj,:))*2/(xmax-xmin);
                            cdRdy(t) = PHTelem{patchIndex}(i).modifiedC(t,:) * squeeze(dBdvDE(ii,jj,:))*2/(ymax-ymin);
                        end
                    end
                    
                    coord{jj,ii} = [coord_pt(1), coord_pt(2)];
                    stress = ghole(coord_pt(1), coord_pt(2), rad, tx);
                    
                    %                    Solve for first derivatives in global coordinates
                    if abs(det(dxdxi))>1e-3
                        dR = dxdxi\[cdRdx';cdRdy'];
                    else
                        dR = pinv(dxdxi)*[cdRdx';cdRdy'];
                    end
                    
                    B = zeros(2*nument,3);
                    B(1:2:2*nument-1,1) = dR(1,:);
                    B(2:2:2*nument,2) = dR(2,:);
                    B(1:2:2*nument-1,3) = dR(2,:);
                    B(2:2:2*nument,3) = dR(1,:);
                    
                    %calculate displacement values
                    dispmatx(jj,ii) = cR'*sol0(scrt_x);
                    dispmaty(jj,ii) = cR'*sol0(scrt_y);
                    
                    %calculate the error in stress values
                    stressvect{jj,ii} = stress - Cmat*B'*sol0(dscrtx);
                    
                    
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
            
            sigmacoord((elementCounter-1)*4+1, :) = stressvect{1,1}';
            sigmacoord((elementCounter-1)*4+2, :) = stressvect{1,2}';
            sigmacoord((elementCounter-1)*4+3, :) = stressvect{2,2}';
            sigmacoord((elementCounter-1)*4+4, :) = stressvect{2,1}';
            
            hold on
        end
    end
end
factor = 5000;

figure
%set(gca,'fontsize',40)
%get(gca)
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sigmacoord(:,1), 'EdgeColor','none', 'facecolor','interp')
view(0,90)
colorbar
title('Displacements and error in \sigma_{11}')
axis off
drawnow

figure
%set(gca,'fontsize',40)
%get(gca)
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sigmacoord(:,3), 'EdgeColor','none', 'facecolor','interp')
view(0,90)
colorbar
title('Displacements and error in \sigma_{12}')
axis off
drawnow

figure
%set(gca,'fontsize',40)
%get(gca)
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sigmacoord(:,2), 'EdgeColor','none', 'facecolor','interp')
view(0,90)
colorbar
%title('Displacements and error in \sigma_{22}')
axis off
drawnow


