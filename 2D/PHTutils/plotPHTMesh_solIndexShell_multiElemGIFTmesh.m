function  plotPHTMesh_solIndexShell_multiElemGIFTmesh( PHUTelem,GIFTmesh,p)
%plot solIndex of selected patch
%patch A - blue
% patch B -red
colorArray = {'red', 'magenta', 'blue', 'green'};

numPts = p+1; %number of plot points to use on each edge
space=0.3;
uref = linspace(-1+space,1-space,numPts);
vref = linspace(-1+space,1-space,numPts);
urefEdge = linspace(-1,1,numPts);
vrefEdge = linspace(-1,1,numPts);


for patchIndex=1:length(PHUTelem)
    for elemIndex=1:length(PHUTelem{patchIndex})
        colorIndex = rem((patchIndex-1),4)+1;
        if isempty(PHUTelem{patchIndex}(elemIndex).children)
            
            
            pxmin = PHUTelem{patchIndex}(elemIndex).vertex(1);
            pymin = PHUTelem{patchIndex}(elemIndex).vertex(2);
            pxmax = PHUTelem{patchIndex}(elemIndex).vertex(3);
            pymax = PHUTelem{patchIndex}(elemIndex).vertex(4);
            
            count=0;
            solIndex=PHUTelem{patchIndex}(elemIndex).solIndex;
            for jj=1:numPts
                for ii=1:numPts
                    count=count+1;
                    
                    
                    [coord, ~] = paramMapShell_multiElemGIFTmesh(elemIndex,GIFTmesh{patchIndex},  uref(ii),vref(jj), pxmin, pymin, pxmax, pymax);
                    %coord=paramMap( GIFTmesh{patchIndex},uref(ii),vref(jj), xmin, ymin, xmax, ymax);
                    if solIndex(count)~=0
                        text(coord(1), coord(2),coord(3),num2str(solIndex(count)), 'Color', colorArray{colorIndex},'FontSize', 8)
                        
                    end
                end
            end
            
            coord_south = zeros(numPts,3);
            coord_north = zeros(numPts,3);
            coord_east = zeros(numPts,3);
            coord_west = zeros(numPts,3);
            
            for ii=1:numPts
                
                south =paramMapShell_multiElemGIFTmesh(elemIndex, GIFTmesh{patchIndex},urefEdge(ii),vrefEdge(1), pxmin, pymin, pxmax, pymax);
                coord_south(ii,:) =[south(1),south(2),south(3)];
                
                north=paramMapShell_multiElemGIFTmesh( elemIndex,GIFTmesh{patchIndex},urefEdge(ii),vrefEdge(end), pxmin, pymin, pxmax, pymax);
                coord_north(ii,:) =[north(1),north(2),north(3)];
                
                east= paramMapShell_multiElemGIFTmesh( elemIndex,GIFTmesh{patchIndex},urefEdge(end),vrefEdge(ii), pxmin, pymin, pxmax, pymax);
                coord_east(ii,:) =[east(1),east(2),east(3)];
                
                west= paramMapShell_multiElemGIFTmesh(elemIndex, GIFTmesh{patchIndex},urefEdge(1),vrefEdge(ii), pxmin, pymin, pxmax, pymax);
                coord_west(ii,:) =[west(1),west(2),west(3)];
                
            end
            
            line(coord_south(:,1),coord_south(:,2),coord_south(:,3),'Color', colorArray{colorIndex})
            line(coord_east(:,1),coord_east(:,2),coord_east(:,3),'Color', colorArray{colorIndex})
            line(coord_north(:,1),coord_north(:,2),coord_north(:,3),'Color', colorArray{colorIndex})
            line(coord_west(:,1),coord_west(:,2),coord_west(:,3),'Color', colorArray{colorIndex})
            hold on
            axis equal
            %drawnow
            
        end
    end
end






