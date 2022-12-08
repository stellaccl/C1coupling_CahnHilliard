function [ stiff, rhs, bcdof, bcval] = imposeDirichletBiharmonic_annulus(stiff, rhs, PHTelem,GIFTmesh, p, q ,type2Basis)

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

[nodesDE ] = getNodes(p+1,q+1);

right_nodesDE=nodesDE.right_nodes1;
right_nodes2DE=nodesDE.right_nodes2;

left_nodesDE=nodesDE.left_nodes1;
left_nodes2DE=nodesDE.left_nodes2;

down_nodesDE=nodesDE.down_nodes1;
down_nodes2DE=nodesDE.down_nodes2;

up_nodesDE=nodesDE.up_nodes1;
up_nodes2DE=nodesDE.up_nodes2;

%set the boundary degree of freedom and elements from the 1st patch
for indexPatch = 1:length(PHTelem)
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            
            numBasis=length(PHTelem{indexPatch}(i).nodesGlobal);
            
            if isempty(PHTelem{indexPatch}(i).neighbor_right)
                 
                nodes=[right_nodes,right_nodes2];
               
                for ii=1:length(nodes)
                    if nodes(ii)<=numBasis
                        bcdof_right = [bcdof_right, PHTelem{indexPatch}(i).nodesGlobal(nodes(ii))];
                    end
                end
%                 bcdof_right = [bcdof_right, PHTelem{indexPatch}(i).nodesGlobal(right_nodes), PHTelem{indexPatch}(i).nodesGlobal(right_nodes2)];
            end
            
            if isempty(PHTelem{indexPatch}(i).neighbor_left)
                
                nodes=[left_nodes,left_nodes2];
                
                for ii=1:length(nodes)
                    if nodes(ii)<=numBasis
                        bcdof_left = [bcdof_left, PHTelem{indexPatch}(i).nodesGlobal(nodes(ii))];
                    end
                end
%                 bcdof_left = [bcdof_left, PHTelem{indexPatch}(i).nodesGlobal(left_nodes), PHTelem{indexPatch}(i).nodesGlobal(left_nodes2)];
            end
            
        end
    end
end


bcdof = unique([bcdof_left, bcdof_right]);
bcdof = setdiff(bcdof,type2Basis);
bcval = zeros(size(bcdof));


[stiff,rhs]=feaplyc2sym(stiff, rhs, bcdof, bcval);
end

