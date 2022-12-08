function  plotPHTMesh_nodesGlobalShell_multiElemGIFTmesh_selectedPatch( PHUTelem,GIFTmesh,p,patch)
%plot solIndex of selected patch
%patch A - blue
% patch B -red
colorArray = {'red', 'magenta', 'blue', 'green'};
%colorArray = {'blue', 'red', 'green', 'cyan', 'magenta', 'black'};
numPts = p+1; %number of plot points to use on each edge
space=0.3;
uref = linspace(-1+space,1-space,numPts);
vref = linspace(-1+space,1-space,numPts);
urefEdge = linspace(-1,1,numPts);
vrefEdge = linspace(-1,1,numPts);


for patchIndex=patch
    for elemIndex=1:length(PHUTelem{patchIndex})
        colorIndex = rem((patchIndex-1),4)+1;
        if isempty(PHUTelem{patchIndex}(elemIndex).children)

            pxmin = PHUTelem{patchIndex}(elemIndex).vertex(1);
            pymin = PHUTelem{patchIndex}(elemIndex).vertex(2);
            pxmax = PHUTelem{patchIndex}(elemIndex).vertex(3);
            pymax = PHUTelem{patchIndex}(elemIndex).vertex(4);
            
            count=0;
            nodesGlobal=PHUTelem{patchIndex}(elemIndex).nodesGlobal;
            for jj=1:numPts
                for ii=1:numPts
                    count=count+1;
                    [ coord, ~] = paramMapShell_multiElemGIFTmesh(elemIndex,GIFTmesh{patchIndex},  uref(ii),vref(jj), pxmin, pymin, pxmax, pymax);
                    text(coord(1), coord(2),coord(3),num2str(nodesGlobal(count)), 'Color', colorArray{colorIndex} ,'FontSize', 8)
                end
            end

            coord_south = zeros(numPts,3);
            coord_east = zeros(numPts,3);
            coord_north = zeros(numPts,3);
            coord_west = zeros(numPts,3);
            
            for j=1:numPts
                coord_south(j,:) = paramMapShell_multiElemGIFTmesh(elemIndex,GIFTmesh{patchIndex}, urefEdge(j), -1, pxmin, pymin, pxmax, pymax );
                coord_east(j,:) = paramMapShell_multiElemGIFTmesh(elemIndex,GIFTmesh{patchIndex}, 1, vrefEdge(j), pxmin, pymin, pxmax, pymax );
                coord_north(j,:) = paramMapShell_multiElemGIFTmesh(elemIndex,GIFTmesh{patchIndex}, -urefEdge(j), 1, pxmin, pymin, pxmax, pymax );
                coord_west(j,:) = paramMapShell_multiElemGIFTmesh(elemIndex,GIFTmesh{patchIndex}, -1, -vrefEdge(j), pxmin, pymin, pxmax, pymax );
            end
            
            coords = [coord_south; coord_east; coord_north; coord_west];
            hold on
            line(coords(:,1), coords(:,2), coords(:,3), 'Color', colorArray{colorIndex})
            
        end
    end
end

hold on
drawnow
axis equal




