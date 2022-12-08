function [ stiff, rhs, bcdof, bcval] = imposeDirichletBiharmonic_squareWithHoleModified_degElev(stiff, rhs, PHTelem, PHTelemDE,GIFTmesh, p, q ,type2Basis)

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


rightNodes=[right_nodes,right_nodes2];
leftNodes=[left_nodes,left_nodes2];
upNodes=[up_nodes,up_nodes2];
downNodes=[down_nodes,down_nodes2];

rightNodesDE=[right_nodesDE,right_nodes2DE];
leftNodesDE=[left_nodesDE,left_nodes2DE];
upNodesDE=[up_nodesDE,up_nodes2DE];
downNodesDE=[down_nodesDE,down_nodes2DE];

indexPatch=1;
for i=1:length(PHTelem{indexPatch})
    
    if isempty(PHTelem{indexPatch}(i).neighbor_up)
        bcdof_up = [bcdof_up, PHTelem{indexPatch}(i).remainNodes(upNodes)];
        bcdof_up = [bcdof_up, PHTelemDE{indexPatch}(i).remainNodes(upNodesDE)];
    end
    
    if isempty(PHTelem{indexPatch}(i).neighbor_down)
        bcdof_down = [bcdof_down, PHTelem{indexPatch}(i).remainNodes(downNodes)];
        bcdof_down = [bcdof_down, PHTelemDE{indexPatch}(i).remainNodes(downNodesDE)];
    end
    
end

indexPatch=2;
for i=1:length(PHTelem{indexPatch})
    
    if isempty(PHTelem{indexPatch}(i).neighbor_up)
        bcdof_up = [bcdof_up, PHTelem{indexPatch}(i).remainNodes(upNodes)];
        bcdof_up = [bcdof_up, PHTelemDE{indexPatch}(i).remainNodes(upNodesDE)];
    end
    
    if isempty(PHTelem{indexPatch}(i).neighbor_left)
        bcdof_left = [bcdof_left, PHTelem{indexPatch}(i).remainNodes(leftNodes)];
        bcdof_left = [bcdof_left, PHTelemDE{indexPatch}(i).remainNodes(leftNodesDE)];
    end
    
end

indexPatch=3;
for i=1:length(PHTelem{indexPatch})
    
    if isempty(PHTelem{indexPatch}(i).neighbor_right)
        bcdof_right = [bcdof_right, PHTelem{indexPatch}(i).remainNodes(rightNodes)];
        bcdof_right = [bcdof_right, PHTelemDE{indexPatch}(i).remainNodes(rightNodesDE)];
    end
    
    if isempty(PHTelem{indexPatch}(i).neighbor_left)
        bcdof_left = [bcdof_left, PHTelem{indexPatch}(i).remainNodes(leftNodes)];
        bcdof_left = [bcdof_left, PHTelemDE{indexPatch}(i).remainNodes(leftNodesDE)];
    end
    
end

indexPatch=4;
for i=1:length(PHTelem{indexPatch})
    
    if isempty(PHTelem{indexPatch}(i).neighbor_down)
        bcdof_down = [bcdof_down, PHTelem{indexPatch}(i).remainNodes(downNodes)];
        bcdof_down = [bcdof_down, PHTelemDE{indexPatch}(i).remainNodes(downNodesDE)];
    end
    
    if isempty(PHTelem{indexPatch}(i).neighbor_left)
        bcdof_left = [bcdof_left, PHTelem{indexPatch}(i).remainNodes(leftNodes)];
        bcdof_left = [bcdof_left, PHTelemDE{indexPatch}(i).remainNodes(leftNodesDE)];
    end
    
end


indexPatch=5;
for i=1:length(PHTelem{indexPatch})
    
    if isempty(PHTelem{indexPatch}(i).neighbor_down)
        bcdof_down = [bcdof_down, PHTelem{indexPatch}(i).remainNodes(downNodes)];
        bcdof_down = [bcdof_down, PHTelemDE{indexPatch}(i).remainNodes(downNodesDE)];
        
    end
    
    if isempty(PHTelem{indexPatch}(i).neighbor_up)
        bcdof_up = [bcdof_up, PHTelem{indexPatch}(i).remainNodes(upNodes)];
        bcdof_up = [bcdof_up, PHTelemDE{indexPatch}(i).remainNodes(upNodesDE)];
    end
    
end

indexPatch=6;
for i=1:length(PHTelem{indexPatch})
    
    if isempty(PHTelem{indexPatch}(i).neighbor_down)
        bcdof_down = [bcdof_down, PHTelem{indexPatch}(i).remainNodes(downNodes)];
        bcdof_down = [bcdof_down, PHTelemDE{indexPatch}(i).remainNodes(downNodesDE)];
    end
    
    if isempty(PHTelem{indexPatch}(i).neighbor_right)
        bcdof_right = [bcdof_right, PHTelem{indexPatch}(i).remainNodes(rightNodes)];
        bcdof_right = [bcdof_right, PHTelemDE{indexPatch}(i).remainNodes(rightNodesDE)];
    end
    
end


indexPatch=7;
for i=1:length(PHTelem{indexPatch})
    
    if isempty(PHTelem{indexPatch}(i).neighbor_left)
        bcdof_left = [bcdof_left, PHTelem{indexPatch}(i).remainNodes(leftNodes)];
        bcdof_left = [bcdof_left, PHTelemDE{indexPatch}(i).remainNodes(leftNodesDE)];
    end
    
    if isempty(PHTelem{indexPatch}(i).neighbor_right)
        bcdof_right = [bcdof_right, PHTelem{indexPatch}(i).remainNodes(rightNodes)];
        bcdof_right = [bcdof_right, PHTelemDE{indexPatch}(i).remainNodes(rightNodesDE)];
    end
    
end

indexPatch=8;
for i=1:length(PHTelem{indexPatch})
    
    if isempty(PHTelem{indexPatch}(i).neighbor_up)
        bcdof_up = [bcdof_up, PHTelem{indexPatch}(i).remainNodes(upNodes)];
        bcdof_up = [bcdof_up, PHTelemDE{indexPatch}(i).remainNodes(upNodesDE)];
    end
    
    if isempty(PHTelem{indexPatch}(i).neighbor_right)
        bcdof_right = [bcdof_right, PHTelem{indexPatch}(i).remainNodes(rightNodes)];
        bcdof_right = [bcdof_right, PHTelemDE{indexPatch}(i).remainNodes(rightNodesDE)];
    end
    
end


bcdof = unique([bcdof_up,bcdof_down,bcdof_right,bcdof_left]);
bcdof=nonzeros(bcdof);
bcdof = setdiff(bcdof,type2Basis);

bcdof =bcdof';
bcval = zeros(size(bcdof));


[stiff,rhs]=feaplyc2sym(stiff, rhs, bcdof, bcval);
end

