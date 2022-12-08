function [nodes ] = getNodes(p,q)
nodes=struct;

nodes.right_nodes1 = (p+1):(p+1):(p+1)*(q+1);
nodes.right_nodes2 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
nodes.right_nodes3 = (p+1)-2:(p+1):(p+1)*(q+1)-2;

nodes.left_nodes1 = 1:(p+1):(1+(p+1)*q);
nodes.left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;
nodes.left_nodes3 = 3:(p+1):(1+(p+1)*q)+2;

nodes.down_nodes1=1:(p+1);
nodes.down_nodes2=nodes.down_nodes1+(p+1);
nodes.down_nodes3=nodes.down_nodes2+(p+1);

nodes.up_nodes1=(p+1)*(q+1)-p:(p+1)*(q+1);
nodes.up_nodes2=nodes.up_nodes1-(p+1);
nodes.up_nodes3=nodes.up_nodes2-(p+1);
end

