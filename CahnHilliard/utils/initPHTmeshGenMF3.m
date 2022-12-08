function [ PHTelem, dimBasis, tempDimBasis] = initPHTmeshGenMF3( p,q, numElemU, numElemV )
%initialize the PHT geometry on coarse mesh, with numElemU x numElemV elements
%use non-open knot vector

alpha = floor((p-1)/2);
beta = floor((q-1)/2);

knotU = [zeros(1,p+1), (1:numElemU)/numElemU, ones(1,p)];
knotV = [zeros(1,q+1), (1:numElemV)/numElemV, ones(1,q)];

%repeat the interior knots p-alpha times
rep_knotU = linspace(0,1,numElemU+1);
rep_knotU = rep_knotU(2:end-1);
rep_knotU = repmat(rep_knotU,1,p-alpha-1);

rep_knotV = linspace(0,1,numElemV+1);
rep_knotV = rep_knotV(2:end-1);
rep_knotV = repmat(rep_knotV,1,q-beta-1);

knotU = sort([knotU, rep_knotU]);
knotV = sort([knotV, rep_knotV]);

tempVector=linspace(0,1,numElemU+1);

diff=tempVector(2)-tempVector(1);
tempVector = repelem(tempVector,2);
knotU_temp=zeros(1,size(tempVector,2)+(p+1)*2);
knotU_temp(1:p+1)=tempVector(1)-diff;
knotU_temp(end-p:end)=tempVector(end)+diff;
knotU_temp(p+2:end-(p+1))=tempVector;

tempVector=linspace(0,1,numElemV+1);

diff=tempVector(2)-tempVector(1);
tempVector = repelem(tempVector,2);
knotV_temp=zeros(1,size(tempVector,2)+(p+1)*2);
knotV_temp(1:p+1)=tempVector(1)-diff;
knotV_temp(end-p:end)=tempVector(end)+diff;
knotV_temp(p+2:end-(p+1))=tempVector;

[C_u_temp, ~] = bezierExtraction(knotU_temp,p);
[C_v_temp, ~] = bezierExtraction(knotV_temp,q);

C_u = C_u_temp(:,:,2:end-1);
C_v = C_v_temp(:,:,2:end-1);

%compute the number of basis functions in each direction
lenU = length(knotU)-p-1;
lenV = length(knotV)-q-1;


%tempDimBasis =one element less on u and v directions.
tempLenU=lenU-2;
tempLenV=lenV-2;
tempDimBasis=tempLenU*tempLenV;

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
        if (knotU(i+1)>knotU(i)) && (knotV(j+1) >knotV(j)) %&& (knotU(i)>=0) && (knotU(i+1)<=1) && (knotV(i)>=0) && (knotV(i+1)<=1)  %the knotspan has non-zero area
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
            PHTelem(elementCounter).nodesGlobal = PHTelem(elementCounter).nodes;
            PHTelem(elementCounter).level = 0;
            PHTelem(elementCounter).isActive = 1;
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
        PHTelem(elementIndex).modifiedC = PHTelem(elementIndex).C;
        
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




