function [ PHTelem, dimBasis] = initPHTmesh_quadratic( p,q ,numElemU,numElemV)


knotU = [zeros(1,p+1), (1:numElemU)/numElemU, ones(1,p)];
knotV = [zeros(1,q+1), (1:numElemV)/numElemV, ones(1,q)];


[C_u, ~] = bezierExtraction(knotU,p);
[C_v, ~] = bezierExtraction(knotV,q);

lenU = length(knotU)-p-1;
lenV = length(knotV)-q-1;

dimBasis = lenU*lenV;
numElements = numElemU*numElemV;

PHTelem = struct;
%initialize the neighbor connectivity lists
PHTelem.neighbor_left = [];
PHTelem.neighbor_right = [];
PHTelem.neighbor_down = [];
PHTelem.neighbor_up = [];

%loop through each element and compute the element-node connectivities
nument = (p+1)*(q+1);
elementCounter = 0;

for j=1:length(knotV)-1
    for i=1:length(knotU)-1
        if (knotU(i+1)>knotU(i)) && (knotV(j+1) >knotV(j))  %the knotspan has non-zero area
            elementCounter = elementCounter + 1;
            PHTelem(elementCounter).parent = [];
            PHTelem(elementCounter).children = [];
            PHTelem(elementCounter).vertex = [knotU(i), knotV(j), knotU(i+1), knotV(j+1)];
            
            tcount = 0;
            currow = zeros(1, nument);
            %now we add the nodes from i-p...i in the u
            %direction, j-q...j in the v direction
            
            for t2=j-q:j
                for t1 = i-p:i
                    tcount = tcount + 1;
                    currow(tcount) = t1+(t2-1)*lenU;
                end
            end
            
            PHTelem(elementCounter).nodes=currow;
            PHTelem(elementCounter).level = 0;
            
        end
    end
end


%loop through each element and compute the neighbor lists and Bezier
%extraction operators
indexMatrix = permute(reshape(1:numElements, numElemU, numElemV),[2,1]);

for j=1:numElemV
    for i=1:numElemU
        elementIndex = indexMatrix(j,i);
        PHTelem(elementIndex).C = kron(C_v(:,:,j),C_u(:,:,i));
        
        if i>1
            PHTelem(elementIndex).neighbor_left = indexMatrix(j,i-1);
        end
        if i<numElemU
            PHTelem(elementIndex).neighbor_right = indexMatrix(j,i+1);
        end
        
        if j>1
            PHTelem(elementIndex).neighbor_down = indexMatrix(j-1,i);
        end
        if j<numElemV
            PHTelem(elementIndex).neighbor_up = indexMatrix(j+1,i);
        end
    end
end

end

