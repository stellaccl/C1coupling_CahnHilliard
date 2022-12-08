function [ PHUTelem,dimBasis] = reassignNodesPeriodicBoundary( PHUTelem,numPatches,p,q,sizeBasis )
%change basis nodes
%after zipConforming

right_nodes1 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
right_nodes2 = (p+1):(p+1):(p+1)*(q+1);

left_nodes1 = 1:(p+1):(1+(p+1)*q);
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

down_nodes1= 1:p+1;
down_nodes2=down_nodes1+(p+1);

up_nodes2=(p+1)*(q+1)-p:(p+1)*(q+1);
up_nodes1=up_nodes2-(p+1);

%find  edge elem (up down left right)
up_edge=[];
down_edge=[];
left_edge=[];
right_edge=[];

for indexPatch=1:numPatches
    for indexElem=1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(indexElem).children)
            if isempty(PHUTelem{indexPatch}(indexElem).neighbor_left) && (indexPatch==1)
                left_edge=[left_edge;indexPatch,indexElem];
            end
            
            if isempty(PHUTelem{indexPatch}(indexElem).neighbor_right) && (indexPatch==2)
                right_edge=[right_edge;indexPatch,indexElem];
            end
            
            if isempty(PHUTelem{indexPatch}(indexElem).neighbor_up)
                up_edge=[up_edge;indexPatch,indexElem];
            end
            
            if isempty(PHUTelem{indexPatch}(indexElem).neighbor_down)
                down_edge=[down_edge;indexPatch,indexElem];
            end
        end
    end
    
end

removedNodes=[];
for i=1:size(up_edge,1)
    
    removedNodes=[removedNodes,PHUTelem{up_edge(i,1)}(up_edge(i,2)).nodesGlobal(up_nodes1)];
    PHUTelem{up_edge(i,1)}(up_edge(i,2)).nodesGlobal(up_nodes1)=PHUTelem{down_edge(i,1)}(down_edge(i,2)).nodesGlobal(down_nodes1);
 %   PHUTelem{up_edge(i,1)}(up_edge(i,2)).C(up_nodes1,:)=PHUTelem{down_edge(i,1)}(down_edge(i,2)).C(down_nodes1,:);
     
    removedNodes=[removedNodes,PHUTelem{up_edge(i,1)}(up_edge(i,2)).nodesGlobal(up_nodes2)];
    PHUTelem{up_edge(i,1)}(up_edge(i,2)).nodesGlobal(up_nodes2)=PHUTelem{down_edge(i,1)}(down_edge(i,2)).nodesGlobal(down_nodes2);
 %   PHUTelem{up_edge(i,1)}(up_edge(i,2)).C(up_nodes2,:)=PHUTelem{down_edge(i,1)}(down_edge(i,2)).C(down_nodes2,:);
end


for i=1:size(right_edge,1)
    
    removedNodes=[removedNodes,PHUTelem{right_edge(i,1)}(right_edge(i,2)).nodesGlobal(right_nodes1)];
    PHUTelem{right_edge(i,1)}(right_edge(i,2)).nodesGlobal(right_nodes1)=PHUTelem{left_edge(i,1)}(left_edge(i,2)).nodesGlobal(left_nodes1);
  %  PHUTelem{right_edge(i,1)}(right_edge(i,2)).C(right_nodes1,:)=PHUTelem{left_edge(i,1)}(left_edge(i,2)).C(left_nodes1,:);
    
    removedNodes=[removedNodes,PHUTelem{right_edge(i,1)}(right_edge(i,2)).nodesGlobal(right_nodes2)];
    PHUTelem{right_edge(i,1)}(right_edge(i,2)).nodesGlobal(right_nodes2)=PHUTelem{left_edge(i,1)}(left_edge(i,2)).nodesGlobal(left_nodes2);
 %   PHUTelem{right_edge(i,1)}(right_edge(i,2)).C(right_nodes2,:)=PHUTelem{left_edge(i,1)}(left_edge(i,2)).C(left_nodes2,:);
end

removedNodes=unique(removedNodes);


%re-assign global nodes
temp_nodes=1:sizeBasis;
replacePattern = temp_nodes;
temp_nodes(removedNodes)=[];
dimBasis=length(temp_nodes);
replacePattern(temp_nodes) = 1:dimBasis;


for indexPatch=1:numPatches
    for indexElem=1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(indexElem).children)
            PHUTelem{indexPatch}(indexElem).nodesGlobal = replacePattern(PHUTelem{indexPatch}(indexElem).nodesGlobal);
        end
    end
end



