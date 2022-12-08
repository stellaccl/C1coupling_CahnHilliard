function  testPlotBasisPhys_GIFTmesh( PHUTelem,GIFTmesh,p,q,type2Basis)
%plots the NEW (type 2) basis function in the PHTelem structure in the physical space
%uses the same color for all elements in a patch to allow plotting with
%more elements


colorArray = {'blue', 'red', 'green', 'cyan', 'magenta', 'black'};
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

%nument = length(PHUTelem{indexPatch}(indexElem).nodes);
%cpts
cR = zeros(numPts, numPts);
c_dRdu = zeros(numPts, numPts);
c_dRdv = zeros(numPts, numPts);

coordX = zeros(numPts, numPts);
coordY = zeros(numPts, numPts);

figure('position',[10 10 800 800]);

numPatches = length(PHUTelem);
for indexBasis=type2Basis
    %for indexBasis=dimBasis-numType2Basis+1:dimBasis
    indexBasis
    
    for indexPatch = 1:numPatches
       % colorIndex = rem((indexPatch-1),6)+1;
       
       for indexElem = 1:length(PHUTelem{indexPatch})
            if isempty(PHUTelem{indexPatch}(indexElem).children)
                if ismember(indexBasis,PHUTelem{indexPatch}(indexElem).nodesGlobal)
                colorIndex = rem((indexElem-1),6)+1;
             
                    localIndex = (PHUTelem{indexPatch}(indexElem).nodesGlobal==indexBasis);
                    
                    %PHUTelem{indexPatch}(indexElem).vertex
                    xmin = PHUTelem{indexPatch}(indexElem).vertex(1);
                    xmax = PHUTelem{indexPatch}(indexElem).vertex(3);
                    ymin = PHUTelem{indexPatch}(indexElem).vertex(2);
                    ymax = PHUTelem{indexPatch}(indexElem).vertex(4);
                    
                    for jj=1:numPts
                        for ii=1:numPts
                            
                            
                            cR(jj,ii) = (PHUTelem{indexPatch}(indexElem).modifiedC(localIndex,:))*squeeze(R(ii,jj,:));
                            [ coord, dxdxi] = paramMap( GIFTmesh{indexPatch}, uref(ii), vref(jj), xmin, ymin, xmax, ymax);
                            
                            % dxdxi=dxdxi';
                            
                            dFloc(1,:) = (PHUTelem{indexPatch}(indexElem).modifiedC(localIndex,:))*squeeze(dRdu(ii,jj,:));
                            dFloc(2,:) = (PHUTelem{indexPatch}(indexElem).modifiedC(localIndex,:))*squeeze(dRdv(ii,jj,:));
                            
                            dFglob = dxdxi\dFloc;
                            
                            c_dRdu(jj,ii) = dFglob(1);
                            c_dRdv(jj,ii) = dFglob(2);
                            
                            coordX(jj,ii) = coord(1);
                            coordY(jj,ii) = coord(2);
                        end
                    end
                    
                    %
                    %                     if indexPatch ==1
                    %                         color='green';
                    %                     else
                    %                         color='blue';
                    %                     end
                    %
                    %
                    
                    subplot(3,1,1)
                    surf(coordX,coordY,cR,'FaceColor',colorArray{colorIndex})%,'EdgeColor','none','FaceLighting','phong')
                    title(['basis functions :',num2str(indexBasis)])
                    drawnow
                    hold on
                    subplot(3,1,2)
                    surf(coordX,coordY,c_dRdu,'FaceColor',colorArray{colorIndex})%,'EdgeColor','none','FaceLighting','phong')
                    title('x-derivatives')
                    drawnow
                    hold on
                    subplot(3,1,3)
                    surf(coordX,coordY,c_dRdv,'FaceColor',colorArray{colorIndex})%,'EdgeColor','none','FaceLighting','phong')
                    title('y-derivatives')
                    drawnow
                    hold on
                    
                end
            end
        end
    end
    
    pause
    close all
    figure('position',[10 10 800 800]);
    
end

camlight
close all
