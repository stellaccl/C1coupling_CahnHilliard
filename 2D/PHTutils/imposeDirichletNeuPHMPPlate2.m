function [ stiff, rhs, bcdof, bcval ] = imposeDirichletNeuPHMPPlate2(stiff, rhs, PHTelem, p, q)
%impose Dirichlet boundary conditions for the simply supported plate problem
%fixed (homogeneous) boundary conditions on the left, right, top, bottom sides
%for one patch only
% 2 nodes at each boundary

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
for i=1:length(PHTelem{1})
    if isempty(PHTelem{1}(i).children)                
        if isempty(PHTelem{1}(i).neighbor_right)
            bcdof_right = [bcdof_right, PHTelem{1}(i).nodesGlobal(right_nodes),PHTelem{1}(i).nodesGlobal(right_nodes2)];            
        end      
        if isempty(PHTelem{1}(i).neighbor_down)
            bcdof_down = [bcdof_down, PHTelem{1}(i).nodesGlobal(down_nodes), PHTelem{1}(i).nodesGlobal(down_nodes2)];            
        end      
        if isempty(PHTelem{1}(i).neighbor_up)
            bcdof_up = [bcdof_up, PHTelem{1}(i).nodesGlobal(up_nodes),PHTelem{1}(i).nodesGlobal(up_nodes2)];            
        end      
        if isempty(PHTelem{1}(i).neighbor_left)
            bcdof_left = [bcdof_left, PHTelem{1}(i).nodesGlobal(left_nodes),PHTelem{1}(i).nodesGlobal(left_nodes2)];            
        end                
    end
end


%remove duplicated entries
bcdof = unique([bcdof_down, bcdof_up, bcdof_left, bcdof_right]);
bcval = zeros(size(bcdof));

%update the stiffness matrix and rhs with Dirichlet boundary values
%bcdof
%bcval
[stiff,rhs]=feaplyc2sym(stiff, rhs, bcdof, bcval);
