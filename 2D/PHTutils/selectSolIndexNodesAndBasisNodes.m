function [nodesA,midNodesA,nodesB,midNodesB] = selectSolIndexNodesAndBasisNodes( edgeA,edgeB,p,q )
%select relevant nodes based on edge
right_nodes2 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
right_nodes1 = (p+1):(p+1):(p+1)*(q+1);

left_nodes1 = 1:(p+1):(1+(p+1)*q);
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

down_nodes1=1:(p+1);
down_nodes2=down_nodes1+(p+1);

up_nodes1=(p+1)*(q+1)-p:(p+1)*(q+1);
up_nodes2=up_nodes1-(p+1);


switch edgeA
    case 1
        midNodesA=down_nodes1;
        nodesA=down_nodes2;
    case 2
        midNodesA=right_nodes1;
        nodesA=right_nodes2;
    case 3
        midNodesA=up_nodes1;
        nodesA=up_nodes2;
    case 4
        midNodesA=left_nodes1;
        nodesA=left_nodes2;
end

switch edgeB
    case 1
        midNodesB=down_nodes1;
        nodesB=down_nodes2;
    case 2
        midNodesB=right_nodes1;
        nodesB=right_nodes2;
    case 3
        midNodesB=up_nodes1;
        nodesB=up_nodes2;
    case 4
        midNodesB=left_nodes1;
        nodesB=left_nodes2;
end

end

