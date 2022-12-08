function [ stiff, rhs, bcdof, bcval] = imposeDirichletBiharmonic_c1(stiff, rhs, PHTelem, p, q ,type2Basis)

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
for indexPatch = 1:length(PHTelem)
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            if isempty(PHTelem{indexPatch}(i).neighbor_right) && indexPatch==2
                
                for ii=1:length(right_nodes)
                    if right_nodes(ii)<=length(PHTelem{indexPatch}(i).nodesGlobal)
                        bcdof_right = [bcdof_right,PHTelem{indexPatch}(i).nodesGlobal(right_nodes(ii))];
                    end
                end
                
                for ii=1:length(right_nodes2)
                    if right_nodes2(ii)<=length(PHTelem{indexPatch}(i).nodesGlobal)
                        bcdof_right = [bcdof_right,PHTelem{indexPatch}(i).nodesGlobal(right_nodes2(ii))];
                    end
                end
                
                
                %bcdof_right = [bcdof_right, PHTelem{indexPatch}(i).nodesGlobal(right_nodes),PHTelem{indexPatch}(i).nodesGlobal(right_nodes2)];
                
            end
            if isempty(PHTelem{indexPatch}(i).neighbor_down)
                
                for ii=1:length(down_nodes)
                    if down_nodes(ii)<=length(PHTelem{indexPatch}(i).nodesGlobal)
                        bcdof_down = [bcdof_down,PHTelem{indexPatch}(i).nodesGlobal(down_nodes(ii))];
                    end
                end
                
                for ii=1:length(down_nodes2)
                    if down_nodes2(ii)<=length(PHTelem{indexPatch}(i).nodesGlobal)
                        bcdof_down = [bcdof_down,PHTelem{indexPatch}(i).nodesGlobal(down_nodes2(ii))];
                    end
                end
                
                %   bcdof_down = [bcdof_down, PHTelem{indexPatch}(i).nodesGlobal(down_nodes), PHTelem{indexPatch}(i).nodesGlobal(down_nodes2)];
                
            end
            if isempty(PHTelem{indexPatch}(i).neighbor_up)
                
                for ii=1:length(up_nodes)
                    if up_nodes(ii)<=length(PHTelem{indexPatch}(i).nodesGlobal)
                        bcdof_up = [bcdof_up,PHTelem{indexPatch}(i).nodesGlobal(up_nodes(ii))];
                    end
                end
                
                for ii=1:length(up_nodes2)
                    if up_nodes2(ii)<=length(PHTelem{indexPatch}(i).nodesGlobal)
                        bcdof_up = [bcdof_up,PHTelem{indexPatch}(i).nodesGlobal(up_nodes2(ii))];
                    end
                end
                %bcdof_up = [bcdof_up, PHTelem{indexPatch}(i).nodesGlobal(up_nodes),PHTelem{indexPatch}(i).nodesGlobal(up_nodes2)];
                
            end
            if isempty(PHTelem{indexPatch}(i).neighbor_left) && indexPatch==1
                
                for ii=1:length(left_nodes)
                    if left_nodes(ii)<=length(PHTelem{indexPatch}(i).nodesGlobal)
                        bcdof_left = [bcdof_left,PHTelem{indexPatch}(i).nodesGlobal(left_nodes(ii))];
                    end
                end
                
                
                 for ii=1:length(left_nodes2)
                    if left_nodes2(ii)<=length(PHTelem{indexPatch}(i).nodesGlobal)
                        bcdof_left = [bcdof_left,PHTelem{indexPatch}(i).nodesGlobal(left_nodes2(ii))];
                    end
                end
                
                
                %bcdof_left = [bcdof_left, PHTelem{indexPatch}(i).nodesGlobal(left_nodes),PHTelem{indexPatch}(i).nodesGlobal(left_nodes2)];
                
            end
        end
    end
end

bcdof = unique([bcdof_down, bcdof_up, bcdof_left, bcdof_right]);

bcdof = setdiff(bcdof,type2Basis);

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

