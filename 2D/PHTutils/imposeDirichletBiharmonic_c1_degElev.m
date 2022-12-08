function [ stiff, rhs, bcdof, bcval] = imposeDirichletBiharmonic_c1_degElev(stiff, rhs, PHTelem, PHTelemDE ,p, q ,type2Basis)

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

for indexPatch = 1:length(PHTelem)
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            if isempty(PHTelem{indexPatch}(i).neighbor_right) && indexPatch==2
                
                nodes=[right_nodes,right_nodes2];
                
                tempBCdof=PHTelem{indexPatch}(i).remainNodes(nodes);
                tempBCdof=sort(tempBCdof);
                tempBCdof=nonzeros(tempBCdof);
                
                bcdof_right = [bcdof_right;tempBCdof];
                
                
                %%%%%%%%%% for additional nodes 
                nodes=[right_nodesDE,right_nodes2DE];
                
                tempBCdof=PHTelemDE{indexPatch}(i).remainNodes(nodes);
                tempBCdof=sort(tempBCdof);
                tempBCdof=nonzeros(tempBCdof);
                
                bcdof_right = [bcdof_right;tempBCdof];
                %%%%%%%%%%%%%%%%%%%
            end
            
            if isempty(PHTelem{indexPatch}(i).neighbor_down)
                
                nodes=[down_nodes,down_nodes2];
                
                tempBCdof=PHTelem{indexPatch}(i).remainNodes(nodes);
                tempBCdof=sort(tempBCdof);
                tempBCdof=nonzeros(tempBCdof);
                
          
                bcdof_down = [bcdof_down;tempBCdof];
                
                %%%%%%%%%% for additional nodes
                nodes=[down_nodesDE,down_nodes2DE];
                
                tempBCdof=PHTelemDE{indexPatch}(i).remainNodes(nodes);
                tempBCdof=sort(tempBCdof);
                tempBCdof=nonzeros(tempBCdof);
                
                bcdof_down = [bcdof_down;tempBCdof];
            end
            
            if isempty(PHTelem{indexPatch}(i).neighbor_up)
                
                nodes=[up_nodes,up_nodes2];
                
                tempBCdof=PHTelem{indexPatch}(i).remainNodes(nodes);
                tempBCdof=sort(tempBCdof);
                tempBCdof=nonzeros(tempBCdof);
                
                bcdof_up = [bcdof_up;tempBCdof];
                
                %%%%%%%%%% for additional nodes
                nodes=[up_nodesDE,up_nodes2DE];
                
                tempBCdof=PHTelemDE{indexPatch}(i).remainNodes(nodes);
                tempBCdof=sort(tempBCdof);
                tempBCdof=nonzeros(tempBCdof);
                
                bcdof_up = [bcdof_up;tempBCdof];

            end
            
            if isempty(PHTelem{indexPatch}(i).neighbor_left) && indexPatch==1
                
                nodes=[left_nodes,left_nodes2];
                
                tempBCdof=PHTelem{indexPatch}(i).remainNodes(nodes);
                tempBCdof=sort(tempBCdof);
                tempBCdof=nonzeros(tempBCdof);
                
                bcdof_left = [bcdof_left;tempBCdof];
                
                 %%%%%%%%%% for additional nodes
                nodes=[left_nodesDE,left_nodes2DE];
                
                tempBCdof=PHTelemDE{indexPatch}(i).remainNodes(nodes);
                tempBCdof=sort(tempBCdof);
                tempBCdof=nonzeros(tempBCdof);
                
                bcdof_left = [bcdof_left;tempBCdof];
            end
        end
    end
end



bcdof = unique([bcdof_down; bcdof_up;bcdof_left; bcdof_right]);

bcdof = setdiff(bcdof,type2Basis);

bcdof=bcdof';
bcval = zeros(size(bcdof));

% temp_bcdof_right=2*bcdof_right-1;
% temp_bcdof_down=2*bcdof_down;
% temp_bcdof_up=2*bcdof_up;
% temp_bcdof_left=2*bcdof_left-1;

%bcdof = [temp_bcdof_right,temp_bcdof_down,temp_bcdof_up,temp_bcdof_left];
% bcdof = unique([bcdof_down, bcdof_up, bcdof_left, bcdof_right]);
% bcval = zeros(size(bcdof));

[stiff,rhs]=feaplyc2sym(stiff, rhs, bcdof, bcval);
end

