function [ stiff, rhs, bcdof, bcval] = imposeDirichletBiharmonic_squareWith2Holes(stiff, rhs, PHTelem,GIFTmesh, p, q ,type2Basis)

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


rightNodes=[right_nodes,right_nodes2];
leftNodes=[left_nodes,left_nodes2];
upNodes=[up_nodes,up_nodes2];
downNodes=[down_nodes,down_nodes2];
%set the boundary degree of freedom and elements from the 1st patch
indexPatch=1;
for i=1:length(PHTelem{indexPatch})
    if isempty(PHTelem{indexPatch}(i).children)
        
        if isempty(PHTelem{indexPatch}(i).neighbor_right)
            bcdof_right = [bcdof_right, PHTelem{indexPatch}(i).remainNodes(rightNodes)];
        end
        
        if isempty(PHTelem{indexPatch}(i).neighbor_left)
            bcdof_left = [bcdof_left, PHTelem{indexPatch}(i).remainNodes(leftNodes)];
        end
        
    end
end


indexPatch=2;
for i=1:length(PHTelem{indexPatch})
    if isempty(PHTelem{indexPatch}(i).children)
        
        if isempty(PHTelem{indexPatch}(i).neighbor_down)
            bcdof_down = [bcdof_down, PHTelem{indexPatch}(i).remainNodes(downNodes)];
        end

    end
end

indexPatch=3;
for i=1:length(PHTelem{indexPatch})
    if isempty(PHTelem{indexPatch}(i).children)
        
        if isempty(PHTelem{indexPatch}(i).neighbor_down)
            bcdof_down = [bcdof_down, PHTelem{indexPatch}(i).remainNodes(downNodes)];
        end
        
        if isempty(PHTelem{indexPatch}(i).neighbor_up)
            bcdof_up = [bcdof_up, PHTelem{indexPatch}(i).remainNodes(upNodes)];
        end
    end
end


indexPatch=4;
for i=1:length(PHTelem{indexPatch})
    if isempty(PHTelem{indexPatch}(i).children)
        
        if isempty(PHTelem{indexPatch}(i).neighbor_down)
            bcdof_down = [bcdof_down, PHTelem{indexPatch}(i).remainNodes(downNodes)];
        end
        
        if isempty(PHTelem{indexPatch}(i).neighbor_left)
            bcdof_left = [bcdof_left, PHTelem{indexPatch}(i).remainNodes(leftNodes)];
        end
    end
end

indexPatch=5;
for i=1:length(PHTelem{indexPatch})
    if isempty(PHTelem{indexPatch}(i).children)
        
        if isempty(PHTelem{indexPatch}(i).neighbor_right)
            bcdof_right = [bcdof_right, PHTelem{indexPatch}(i).remainNodes(rightNodes)];
        end
        
        if isempty(PHTelem{indexPatch}(i).neighbor_left)
            bcdof_left = [bcdof_left, PHTelem{indexPatch}(i).remainNodes(leftNodes)];
        end
    end
end

indexPatch=6;
for i=1:length(PHTelem{indexPatch})
    if isempty(PHTelem{indexPatch}(i).children)
        
        if isempty(PHTelem{indexPatch}(i).neighbor_up)
            bcdof_up = [bcdof_up, PHTelem{indexPatch}(i).remainNodes(upNodes)];
        end
        
        if isempty(PHTelem{indexPatch}(i).neighbor_left)
            bcdof_left = [bcdof_left, PHTelem{indexPatch}(i).remainNodes(leftNodes)];
        end
    end
end



indexPatch=7;
for i=1:length(PHTelem{indexPatch})
    if isempty(PHTelem{indexPatch}(i).children)
        
        if isempty(PHTelem{indexPatch}(i).neighbor_up)
            bcdof_up = [bcdof_up, PHTelem{indexPatch}(i).remainNodes(upNodes)];
        end
        
        if isempty(PHTelem{indexPatch}(i).neighbor_down)
            bcdof_down = [bcdof_down, PHTelem{indexPatch}(i).remainNodes(downNodes)];
        end
    end
end

indexPatch=8;
for i=1:length(PHTelem{indexPatch})
    if isempty(PHTelem{indexPatch}(i).children)
        
        if isempty(PHTelem{indexPatch}(i).neighbor_up)
            bcdof_up = [bcdof_up, PHTelem{indexPatch}(i).remainNodes(upNodes)];
        end

    end
end


indexPatch=9;
for i=1:length(PHTelem{indexPatch})
    if isempty(PHTelem{indexPatch}(i).children)
        
        if isempty(PHTelem{indexPatch}(i).neighbor_up)
            bcdof_up = [bcdof_up, PHTelem{indexPatch}(i).remainNodes(upNodes)];
        end
        
        if isempty(PHTelem{indexPatch}(i).neighbor_down)
            bcdof_down = [bcdof_down, PHTelem{indexPatch}(i).remainNodes(downNodes)];
        end
    end
end

indexPatch=10;
for i=1:length(PHTelem{indexPatch})
    if isempty(PHTelem{indexPatch}(i).children)
        
        if isempty(PHTelem{indexPatch}(i).neighbor_up)
            bcdof_up = [bcdof_up, PHTelem{indexPatch}(i).remainNodes(upNodes)];
        end
        
        if isempty(PHTelem{indexPatch}(i).neighbor_right)
            bcdof_right = [bcdof_right, PHTelem{indexPatch}(i).remainNodes(rightNodes)];
        end
    end
end


indexPatch=11;
for i=1:length(PHTelem{indexPatch})
    if isempty(PHTelem{indexPatch}(i).children)
        
        if isempty(PHTelem{indexPatch}(i).neighbor_left)
            bcdof_left = [bcdof_left, PHTelem{indexPatch}(i).remainNodes(leftNodes)];
        end
        
        if isempty(PHTelem{indexPatch}(i).neighbor_right)
            bcdof_right = [bcdof_right, PHTelem{indexPatch}(i).remainNodes(rightNodes)];
        end
    end
end


indexPatch=12;
for i=1:length(PHTelem{indexPatch})
    if isempty(PHTelem{indexPatch}(i).children)
        
        if isempty(PHTelem{indexPatch}(i).neighbor_down)
            bcdof_down = [bcdof_down, PHTelem{indexPatch}(i).remainNodes(downNodes)];
        end
        
        if isempty(PHTelem{indexPatch}(i).neighbor_right)
            bcdof_right = [bcdof_right, PHTelem{indexPatch}(i).remainNodes(rightNodes)];
        end
    end
end


indexPatch=13;
for i=1:length(PHTelem{indexPatch})
    if isempty(PHTelem{indexPatch}(i).children)
        
        if isempty(PHTelem{indexPatch}(i).neighbor_down)
            bcdof_down = [bcdof_down, PHTelem{indexPatch}(i).remainNodes(downNodes)];
        end
        
        if isempty(PHTelem{indexPatch}(i).neighbor_up)
            bcdof_up = [bcdof_up, PHTelem{indexPatch}(i).remainNodes(upNodes)];
        end
    end
end


bcdof = unique([bcdof_right,bcdof_left,bcdof_up,bcdof_down]);
bcdof = nonzeros(bcdof);
bcdof =bcdof';
bcdof = setdiff(bcdof,type2Basis);
bcval = zeros(size(bcdof));


[stiff,rhs]=feaplyc2sym(stiff, rhs, bcdof, bcval);
end

