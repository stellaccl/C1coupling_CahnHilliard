function  plotPHTShellMeshMP( PHTelem, GIFTmesh )
% plots the PHT mesh for the given shell geometry

numPts = 11; %number of plot points to use on each edge
uref = linspace(-1,1,numPts);
vref = linspace(-1,1,numPts);
%we define colors for up to 6 patches
colorArray = {'red', 'magenta', 'blue', 'green'};

numPatches = length(PHTelem);

for patchIndex = 1:numPatches
    for i=1:length(PHTelem{patchIndex})
         colorIndex = rem((patchIndex-1),4)+1;
        
        if isempty(PHTelem{patchIndex}(i).children)
            
            %store the corners of the element in parameter space for easy
            %access
            pxmin = PHTelem{patchIndex}(i).vertex(1);
            pymin = PHTelem{patchIndex}(i).vertex(2);
            pxmax = PHTelem{patchIndex}(i).vertex(3);
            pymax = PHTelem{patchIndex}(i).vertex(4);
            
            %determine the edges of the element in the physical space
            coord_south = zeros(numPts,3);
            coord_east = zeros(numPts,3);
            coord_north = zeros(numPts,3);
            coord_west = zeros(numPts,3);
            
            for j=1:numPts
               
%                 coord_south(j,:) = paramMapShell_checkPipe(i,GIFTmesh{patchIndex}, uref(j), -1, pxmin, pymin, pxmax, pymax );
%                 coord_east(j,:) = paramMapShell_checkPipe(i,GIFTmesh{patchIndex}, 1, vref(j), pxmin, pymin, pxmax, pymax );
%                 coord_north(j,:) = paramMapShell_checkPipe(i,GIFTmesh{patchIndex}, -uref(j), 1, pxmin, pymin, pxmax, pymax );
%                 coord_west(j,:) = paramMapShell_checkPipe(i,GIFTmesh{patchIndex}, -1, -vref(j), pxmin, pymin, pxmax, pymax );
                
                
                 coord_south(j,:) = paramMapShell(GIFTmesh{patchIndex}, uref(j), -1, pxmin, pymin, pxmax, pymax );
                coord_east(j,:) = paramMapShell(GIFTmesh{patchIndex}, 1, vref(j), pxmin, pymin, pxmax, pymax );
                coord_north(j,:) = paramMapShell(GIFTmesh{patchIndex}, -uref(j), 1, pxmin, pymin, pxmax, pymax );
                coord_west(j,:) = paramMapShell(GIFTmesh{patchIndex}, -1, -vref(j), pxmin, pymin, pxmax, pymax );
%                 
            end
            
            %plot the edges of the element
            
            coords = [coord_south; coord_east; coord_north; coord_west];
            
            line(coords(:,1), coords(:,2), coords(:,3), 'Color', colorArray{ colorIndex})
            
            %write the element number in the middle
            coord_mid = paramMapShell(GIFTmesh{patchIndex}, 0, 0, pxmin, pymin, pxmax, pymax );
            %  coord_mid = paramMapShell_checkPipe(i,GIFTmesh{patchIndex}, 0, 0, pxmin, pymin, pxmax, pymax );
            text(coord_mid(1), coord_mid(2), coord_mid(3), num2str(i), 'Color', colorArray{ colorIndex})
            
            
        end
    end
    drawnow
  
    
end
%axis tight
hold on
drawnow
