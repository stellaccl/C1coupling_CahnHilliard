function  plotPHTMeshMP1D( PHTelem, GIFTmesh )


%plots the elements stored in PHTelem structure array
%supports multipatches

%we define colors for up to 6 patches
colorArray = {'blue', 'red', 'green', 'cyan', 'magenta', 'black'};

numPts = 11; %number of plot points to use on each edge
uref = linspace(-1,1,numPts);

numPatches = length(PHTelem);

for patchIndex = 1:numPatches
    patchIndex
    colorIndex = rem((patchIndex-1),6)+1;
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            
            %store the corners of the element in parameter space for easy
            %access
            pxmin = PHTelem{patchIndex}(i).vertex(1);
            pxmax = PHTelem{patchIndex}(i).vertex(2);
            
            %determine the edges of the element in the physical space
            coords = zeros(numPts,1);
            
            for j=1:numPts
                coords(j) = paramMap1D(GIFTmesh{patchIndex}, uref(j), pxmin, pxmax );
                hold on
                plot(coords(j),0,'.')
            end
            
            %plot the edges of the element
            
           % line(coords(:,1), coords(:,2), 'Color', colorArray{colorIndex}, 'LineWidth',2)
            %line(coords(:,1), coords(:,2), 'Color', 'k', 'LineWidth',2)
            %             %write the element number in the middle
            %              coord_mid = paramMap1D(GIFTmesh{patchIndex}, 0, 0, pxmin, pymin, pxmax, pymax );
            % %             %
            %             text(coord_mid(1), coord_mid(2), num2str(i), 'Color', colorArray{colorIndex})
            %text(coord_mid(1), coord_mid(2), num2str(patchIndex), 'Color', colorArray{colorIndex})
        end
    end
end
axis tight
hold on
drawnow
