function  plotPHTMesh_solIndex( PHUTelem,GIFTmesh,p)
%plot solIndex of selected patch
%patch A - blue
% patch B -red
colorArray = {'blue', 'red', 'green', 'cyan', 'magenta', 'black'};
numPts = p+1; %number of plot points to use on each edge
space=0.3;
uref = linspace(-1+space,1-space,numPts);
vref = linspace(-1+space,1-space,numPts);
urefEdge = linspace(-1,1,numPts);
vrefEdge = linspace(-1,1,numPts);


for patchIndex=1:length(PHUTelem)
     colorIndex = rem((patchIndex-1),6)+1;
    for elemIndex=1:length(PHUTelem{patchIndex})
        if isempty(PHUTelem{patchIndex}(elemIndex).children)
            
            xmin=PHUTelem{patchIndex}(elemIndex).vertex(1);
            ymin=PHUTelem{patchIndex}(elemIndex).vertex(2);
            xmax=PHUTelem{patchIndex}(elemIndex).vertex(3);
            ymax=PHUTelem{patchIndex}(elemIndex).vertex(4);
            
            count=0;
            solIndex=PHUTelem{patchIndex}(elemIndex).solIndex;
            
            for jj=1:numPts
                for ii=1:numPts
                    count=count+1;
                    coord=paramMap( GIFTmesh{patchIndex},uref(ii),vref(jj), xmin, ymin, xmax, ymax);
                    text(coord(1), coord(2),num2str(solIndex(count)), 'Color', colorArray{colorIndex})
                end
            end
            
            coord_south = zeros(numPts,2);
            coord_north = zeros(numPts,2);
            coord_east = zeros(numPts,2);
            coord_west = zeros(numPts,2);
            
            for ii=1:numPts
                
                south = paramMap( GIFTmesh{patchIndex},urefEdge(ii),vrefEdge(1), xmin, ymin, xmax, ymax);
                coord_south(ii,:) =[south(1),south(2)];
                
                north= paramMap( GIFTmesh{patchIndex},urefEdge(ii),vrefEdge(end), xmin, ymin, xmax, ymax);
                coord_north(ii,:) =[north(1),north(2)];
                
                east= paramMap( GIFTmesh{patchIndex},urefEdge(end),vrefEdge(ii), xmin, ymin, xmax, ymax);
                coord_east(ii,:) =[east(1),east(2)];
                
                west= paramMap( GIFTmesh{patchIndex},urefEdge(1),vrefEdge(ii), xmin, ymin, xmax, ymax);
                coord_west(ii,:) =[west(1),west(2)];
                
            end
            
            line(coord_south(:,1),coord_south(:,2),'Color', colorArray{colorIndex}, 'LineWidth',2)
            line(coord_east(:,1),coord_east(:,2),'Color', colorArray{colorIndex}, 'LineWidth',2)
            line(coord_north(:,1),coord_north(:,2),'Color', colorArray{colorIndex}, 'LineWidth',2)
            line(coord_west(:,1),coord_west(:,2),'Color', colorArray{colorIndex}, 'LineWidth',2)
            hold on
            axis equal
            drawnow
            
        end
    end
end






