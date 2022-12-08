function  testPlotBasisPhys_GIFTmesh_flat( PHUTelem,GIFTmesh,p,q,indexBasis)
%plots the NEW (type 2) basis function in the PHTelem structure in the physical space
%uses the same color for all elements in a patch to allow plotting with
%more elements

numPts = 11; %number of plot points to use on each edge
fudge = 0;
uref = linspace(-1+fudge,1-fudge,numPts);
vref = linspace(-1+fudge,1-fudge,numPts);

%1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u,dB_u] = bernstein_basis(uref,p);
[B_v,dB_v] = bernstein_basis(vref,q);

R = zeros(numPts, numPts, (p+1)*(q+1));
dRdu = zeros(numPts, numPts, (p+1)*(q+1));
dRdv = zeros(numPts, numPts, (p+1)*(q+1));
%the derivatives of the 2D Bernstein polynomials at Gauss points on the
%master element
basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        R(:,:,basisCounter) = B_u(:,i)*B_v(:,j)';
        dRdu(:,:,basisCounter) = dB_u(:,i)*B_v(:,j)';
        dRdv(:,:,basisCounter) = B_u(:,i)*dB_v(:,j)';
    end
end
%%%%%%%%%%%%%%%%


%1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u_DE,dB_u_DE] = bernstein_basis(uref,p+1);
[B_v_DE,dB_v_DE] = bernstein_basis(vref,q+1);

R_DE = zeros(numPts, numPts, (p+2)*(q+2));
dRdu_DE = zeros(numPts, numPts, (p+2)*(q+2));
dRdv_DE = zeros(numPts, numPts, (p+2)*(q+2));
%the derivatives of the 2D Bernstein polynomials at Gauss points on the
%master element
basisCounter = 0;
for j=1:q+2
    for i=1:p+2
        basisCounter = basisCounter + 1;
        R_DE(:,:,basisCounter) = B_u_DE(:,i)*B_v_DE(:,j)';
        dRdu_DE(:,:,basisCounter) = dB_u_DE(:,i)*B_v_DE(:,j)';
        dRdv_DE(:,:,basisCounter) = B_u_DE(:,i)*dB_v_DE(:,j)';
    end
end


%%%%%%%%%%%%%%%%%

%nument = length(PHUTelem{indexPatch}(indexElem).nodes);
%cpts
cR = zeros(numPts, numPts);
c_dRdu = zeros(numPts, numPts);
c_dRdv = zeros(numPts, numPts);

coordX = zeros(numPts, numPts);
coordY = zeros(numPts, numPts);

figure('position',[10 10 800 800]);

numPatches = length(PHUTelem);
%for indexBasis=type2Basis
%for indexBasis=dimBasis-numType2Basis+1:dimBasis
%  indexBasis

for indexPatch = 1:numPatches
    
    for indexElem = 1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(indexElem).children)
            if ismember(indexBasis,PHUTelem{indexPatch}(indexElem).nodesGlobal)
                
                localIndex = (PHUTelem{indexPatch}(indexElem).nodesGlobal==indexBasis);
                
                
                %polyDegree=PHUTelem{indexPatch}(indexElem).polyDegree(localIndex);
                
                %PHUTelem{indexPatch}(indexElem).vertex
                xmin = PHUTelem{indexPatch}(indexElem).vertex(1);
                xmax = PHUTelem{indexPatch}(indexElem).vertex(3);
                ymin = PHUTelem{indexPatch}(indexElem).vertex(2);
                ymax = PHUTelem{indexPatch}(indexElem).vertex(4);
                
                for jj=1:numPts
                    for ii=1:numPts
                        
                        %                         if polyDegree==p
                        %
                        cR(jj,ii) = (PHUTelem{indexPatch}(indexElem).modifiedC(localIndex,1:(p+1)*(q+1)))*squeeze(R(ii,jj,:));
                        dFloc(1,:) = (PHUTelem{indexPatch}(indexElem).modifiedC(localIndex,1:(p+1)*(q+1)))*squeeze(dRdu(ii,jj,:));
                        dFloc(2,:) = (PHUTelem{indexPatch}(indexElem).modifiedC(localIndex,1:(p+1)*(q+1)))*squeeze(dRdv(ii,jj,:));
                        
                        %                         elseif polyDegree==p+1
                        %                             cR(jj,ii) = (PHUTelem{indexPatch}(indexElem).modifiedC(localIndex,:))*squeeze(R_DE(ii,jj,:));
                        %                             dFloc(1,:) = (PHUTelem{indexPatch}(indexElem).modifiedC(localIndex,:))*squeeze(dRdu_DE(ii,jj,:));
                        %                             dFloc(2,:) = (PHUTelem{indexPatch}(indexElem).modifiedC(localIndex,:))*squeeze(dRdv_DE(ii,jj,:));
                        %                         end
                        
                        [ coord, dxdxi] = paramMap( GIFTmesh{indexPatch}, uref(ii), vref(jj), xmin, ymin, xmax, ymax);
                        
                        % dxdxi=dxdxi';
                        
                        dFglob = dxdxi\dFloc;
                        
                        c_dRdu(jj,ii) = dFglob(1);
                        c_dRdv(jj,ii) = dFglob(2);
                        
                        coordX(jj,ii) = coord(1);
                        coordY(jj,ii) = coord(2);
                    end
                end
                
                
                if indexPatch ==1
                    color='green';
                else
                    color='blue';
                end
                
                
                
               % subplot(3,1,1)
                 view(0,90)
                surf(coordX,coordY,zeros(size(coordY)),cR,'FaceColor','interp','EdgeColor','none')%,'EdgeColor','none','FaceLighting','phong')
                colorbar
                hold on
                plotPHTMeshMP_boundary( PHUTelem, GIFTmesh )
                title(['basis functions :',num2str(indexBasis)])
                drawnow
                axis equal
                hold on
%                 subplot(3,1,2)
%                  view(0,90)
%                 surf(coordX,coordY,zeros(size(coordY)),c_dRdu,'FaceColor','interp','EdgeColor','none')%,'EdgeColor','none','FaceLighting','phong')
%                 colorbar
%                 hold on
%                 plotPHTMeshMP_boundary( PHUTelem, GIFTmesh )
%                 title('x-derivatives')
%                 drawnow
%                 axis equal
%                 hold on
%                 subplot(3,1,3)
%                  view(0,90)
%                 surf(coordX,coordY,zeros(size(coordY)),c_dRdv,'FaceColor','interp','EdgeColor','none')%,'EdgeColor','none','FaceLighting','phong')
%                 colorbar
%                 hold on
%                 plotPHTMeshMP_boundary( PHUTelem, GIFTmesh )
%                 title('y-derivatives')
%                 drawnow
%                 axis equal
%                 hold on
                
            end
        end
    end
end



% pause
% close all
%figure('position',[10 10 800 800]);

%end
% hold on
% plotPHTMeshMP_boundary( PHUTelem, GIFTmesh )
%camlight
%close all
