function [ stiff, rhs, bcdof, bcval ] = imposeDirichletNeuPHMPCircular(stiff, rhs, PHTelem, p, q, numPatches,type2Basis)
%impose Dirichlet boundary conditions for the simply supported plate problem
%fixed (homogeneous) boundary conditions on the left, right, top, bottom sides
%for one patch only
% 2 nodes at each boundary
% two patches

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
for indexPatch = 1:numPatches
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            
            numBasis=length(PHTelem{indexPatch}(i).nodesGlobal);

            if isempty(PHTelem{indexPatch}(i).neighbor_down) && (indexPatch==1 || indexPatch==2 || indexPatch==3 || indexPatch==4)
                
                nodes=[down_nodes,down_nodes2];
                
                for ii=1:length(nodes)
                    if nodes(ii)<=numBasis
                        bcdof_down = [bcdof_down,PHTelem{indexPatch}(i).nodesGlobal(nodes(ii))];
                    end
                end
                
                % bcdof_down = [bcdof_down, PHTelem{indexPatch}(i).nodesGlobal(down_nodes), PHTelem{indexPatch}(i).nodesGlobal(down_nodes2)];
            end
 
        end
    end
end


%remove duplicated entries
bcdof = unique([bcdof_down, bcdof_up, bcdof_left, bcdof_right]);

%make sure the type2basis are not included in the Dirichlet boundary
%conditions
bcdof = setdiff(bcdof,type2Basis);

bcval = zeros(size(bcdof));

%update the stiffness matrix and rhs with Dirichlet boundary values
%bcdof
%bcval
[stiff,rhs]=feaplyc2sym(stiff, rhs, bcdof, bcval);
