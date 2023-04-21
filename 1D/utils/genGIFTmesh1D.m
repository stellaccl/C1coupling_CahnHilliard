function [ GIFTmesh] = genGIFTmesh1D( knotU, coefs, p, numberElementsU)

dim = 1; %number of physical dimensions

%tolerance for equality tests
toleq = 1e-10;

%the number of control points in the u and v directions
lenU = length(knotU)-p-1;  %number of basis functions in the u direction

numnodes = lenU;
coordinates = zeros(numnodes, dim+1);  %allocate one more column for the weights
for i=1:lenU % for each node in the x direction
    coordinates(i,:) = [coefs(1,i)./coefs(4,i), coefs(4,i)]; %put the (i,j) node in the coordinate array
end
C = zeros((p+1),(p+1),numberElementsU);
[C_u, ~] = bezierExtraction(knotU,p);

elementCounter = 0;
for i=1:numberElementsU
    elementCounter = elementCounter + 1;
    C(:,:,elementCounter) = C_u(:,:,i);
end


elementCounter = 0;
elementVertex = zeros(numberElementsU,2);
elementNode = zeros(numberElementsU, (p+1));


for i=1:length(knotU)-1
    if (abs(knotU(i+1)-knotU(i))>toleq) 
        elementCounter = elementCounter + 1;
        elementVertex(elementCounter, :) = [knotU(i), knotU(i+1)];
%         tcount = 0;
%         currow = zeros(1, (p+1));
        %now we add the nodes from i-p...i in the u direction and
        %j-q...j in the v direction
%         for t2=j-q:j
%             for t1 = i-p:i
%                 tcount = tcount + 1;
  
%             end
%         end
        elementNode(elementCounter,:)=1:2;
    end
end

GIFTmesh.numberElements = numberElementsU;
GIFTmesh.numberElementsU = numberElementsU;


GIFTmesh.p = p;
GIFTmesh.c_net = coordinates;
GIFTmesh.C = C;
GIFTmesh.elementNode = elementNode;
GIFTmesh.elementVertex =elementVertex;

end


