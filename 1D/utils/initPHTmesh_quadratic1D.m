function [ PHTelem, dimBasis] = initPHTmesh_quadratic1D( p,numElemU)

knotU = [zeros(1,p+1), (1:numElemU)/numElemU, ones(1,p)];

[C_u, ~] = bezierExtraction(knotU,p);

lenU = length(knotU)-p-1;


dimBasis = lenU;
numElements = numElemU;

PHTelem = struct;
%initialize the neighbor connectivity lists
PHTelem.neighbor_left = [];
PHTelem.neighbor_right = [];

%loop through each element and compute the element-node connectivities
elementCounter = 0;
initialNodes=1:p+1;
count=0;
for i=1:length(knotU)-1
    if (knotU(i+1)>knotU(i)) %the knotspan has non-zero area
        elementCounter = elementCounter + 1;
        PHTelem(elementCounter).parent = [];
        PHTelem(elementCounter).children = [];
        PHTelem(elementCounter).vertex = [knotU(i), knotU(i+1)];
        PHTelem(elementCounter).nodes=initialNodes+count;
        PHTelem(elementCounter).level = 0;
        count=count+1;
    end
end



%loop through each element and compute the neighbor lists and Bezier
%extraction operators
for i=1:numElemU
    PHTelem(i).C = C_u(:,:,i);
end


end

