function [ stiff, rhs, bcdof, bcval] = imposeDirichletBiharmonic_squareWithHole(stiff, rhs, PHTelem,GIFTmesh, p, q ,type2Basis)

%detect which nodes have support on the left and right boundary
bcdof_left = [];
bcdof_right = [];
bcdof_up = [];
bcdof_down = [];

%define side node indices
down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);

down_nodes2 = p+2:p+(p+2);
right_nodes2 = (p):(p+1):(p+1)*(q+1)-1;
up_nodes2 = (p+1)*q-p:(p+1)*(q+1)-p-1;
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

%set the boundary degree of freedom and elements from the 1st patch
indexPatch=1;
for i=1:length(PHTelem{indexPatch})
    if isempty(PHTelem{indexPatch}(i).children)
        
        if isempty(PHTelem{indexPatch}(i).neighbor_right)
            
            numBasis=length(PHTelem{indexPatch}(i).nodesGlobal);
            nodes=[right_nodes,right_nodes2];
            
            for ii=1:length(nodes)
                if nodes(ii)<=numBasis
                    bcdof_right = [bcdof_right, PHTelem{indexPatch}(i).nodesGlobal(nodes(ii))];
                end
            end
            
        end
        
        if isempty(PHTelem{indexPatch}(i).neighbor_left)
            
            numBasis=length(PHTelem{indexPatch}(i).nodesGlobal);
            nodes=[left_nodes,left_nodes2];
            
            for ii=1:length(nodes)
                if nodes(ii)<=numBasis
                    bcdof_left = [bcdof_left, PHTelem{indexPatch}(i).nodesGlobal(nodes(ii))];
                end
            end
            
        end
        
        
        
    end
end



bcdof = unique([bcdof_right,bcdof_left]);
bcdof = setdiff(bcdof,type2Basis);
bcval = zeros(size(bcdof));


[stiff,rhs]=feaplyc2sym(stiff, rhs, bcdof, bcval);
end

