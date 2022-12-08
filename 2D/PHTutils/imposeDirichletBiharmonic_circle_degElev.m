function [ stiff, rhs, bcdof, bcval] = imposeDirichletBiharmonic_circle_degElev(stiff, rhs,PHTelem,PHTelemDE,GIFTmesh, p, q ,type2Basis,additionalBasis)

%detect which nodes have support on the left and right boundary
bcdof_left = [];
bcdof_right = [];
bcdof_up = [];
bcdof_down = [];

%define side node indices

[nodesDE ] = getNodes(p+1,q+1);

right_nodesDE=nodesDE.right_nodes1;
right_nodes2DE=nodesDE.right_nodes2;

left_nodesDE=nodesDE.left_nodes1;
left_nodes2DE=nodesDE.left_nodes2;

down_nodesDE=nodesDE.down_nodes1;
down_nodes2DE=nodesDE.down_nodes2;

up_nodesDE=nodesDE.up_nodes1;
up_nodes2DE=nodesDE.up_nodes2;

down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);

down_nodes2 = p+2:p+(p+2);
right_nodes2 = (p):(p+1):(p+1)*(q+1)-1;
up_nodes2 = (p+1)*q-p:(p+1)*(q+1)-p-1;
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

%set the boundary degree of freedom and elements from the 1st patch
for indexPatch = 1:4
    for indexElem=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(indexElem).children)
            
            numBasis=length(PHTelem{indexPatch}(indexElem).nodesGlobal);
            
            if isempty(PHTelem{indexPatch}(indexElem).neighbor_down)
                
                nodes=[down_nodes,down_nodes2];
                
                for ii=1:length(nodes)
                    if nodes(ii)<=numBasis
                        bcdof_down = [bcdof_down, PHTelem{indexPatch}(indexElem).nodesGlobal(nodes(ii))];
                    end
                end
            end
            
        end
    end
end

%select nodes from additional basis
for indexPatch = 1:4
    for indexElem=1:length(PHTelemDE{indexPatch})
        
        nodes=[down_nodesDE,down_nodes2DE];
        
        if isempty(PHTelemDE{indexPatch}(indexElem).neighbor_down)
            
            if ~isempty(PHTelemDE{indexPatch}(indexElem).additionalNodes)
                [loc,~]=ismember(PHTelemDE{indexPatch}(indexElem).additionalNodes,nodes);
                bcdof_down = [bcdof_down,PHTelemDE{indexPatch}(indexElem).additionalBasis(loc)];
            end
        end
    end
end


bcdof = unique([bcdof_down]);
bcdof = setdiff(bcdof,type2Basis);
bcval = zeros(size(bcdof));

[stiff,rhs]=feaplyc2sym(stiff, rhs, bcdof, bcval);
end

