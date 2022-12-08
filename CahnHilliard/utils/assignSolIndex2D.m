function  [ PHUTelem ,solIndexCount,patchInfo] = assignSolIndex2D( PHUTelem,patchBoundaries,p,q)
%assign solIndex
%support cross insertion

numBoundaries = size(patchBoundaries,1);
numPatches = length(PHUTelem);

right_nodes2 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
right_nodes1 = (p+1):(p+1):(p+1)*(q+1);

left_nodes1 = 1:(p+1):(1+(p+1)*q);
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

down_nodes1=1:(p+1);
down_nodes2=down_nodes1+(p+1);

up_nodes1=(p+1)*(q+1)-p:(p+1)*(q+1);
up_nodes2=up_nodes1-(p+1);


%create/set nodesGlobal entries in all patches to be equal to local nodes
%entries
for patchIndex = 1:numPatches
    for elemIndex = 1:length(PHUTelem{patchIndex})
        PHUTelem{patchIndex}(elemIndex).solIndex = zeros(1,(p+1)*(q+1));
    end
end



%first row = neighbor patch index
%second row = edge
%third row = neighbor edge
neighborInfo=[];
for boundaryIndex = 1:numBoundaries
    
    patchAList = patchBoundaries{boundaryIndex,1};
    patchB = patchBoundaries{boundaryIndex,2};
    edgeAList = patchBoundaries{boundaryIndex,3};
    edgeBList = patchBoundaries{boundaryIndex,4};
    
    for indexPatch=1:length(patchAList)
        
        patchA = patchAList(indexPatch);
        edgeA = edgeAList(indexPatch);
        edgeB = edgeBList(indexPatch);
        
        neighborInfo=[neighborInfo;patchA,patchB,edgeA,edgeB;patchB,patchA,edgeB,edgeA];
        
    end
    
end
neighborInfo=sortrows(neighborInfo,1);
patchInfo=cell(1,numPatches);
%patch info of each patch contains informaiton of all neighbor patches attached to it
%each row (patchA, patchB,edgeA,edgeB)
for i=1:size(neighborInfo,1)
    patchInfo{neighborInfo(i,1)}=[patchInfo{neighborInfo(i,1)};neighborInfo(i,:)];
end


solIndexCount = 0;

%assign sol index from neighbor patch
for indexPatch=1:numPatches
   neighbor=patchInfo{indexPatch}(:,2);

    for i=1:length(neighbor)
        
        if neighbor(i)<indexPatch
            
            neighborPatch=neighbor(i);
            neighborEdge= patchInfo{indexPatch}(i,4);
            patchEdge= patchInfo{indexPatch}(i,3);
            [~,elemsPatch] = sortEdgeNodesElem( PHUTelem{indexPatch},  patchEdge, p, q );
            [~,elemsNeighbor] = sortEdgeNodesElem( PHUTelem{neighborPatch}, neighborEdge, p, q );
            
            elemsPatch=cell2mat(elemsPatch);
            elemsNeighbor=cell2mat(elemsNeighbor);
            
            switch neighborEdge
                case 1
                    nodesNeighbor=[down_nodes1];
                case 2
                    nodesNeighbor =[right_nodes1];
                case 3
                    nodesNeighbor =[up_nodes1];
                case 4
                    nodesNeighbor =[left_nodes1];
                    
            end
            
            switch patchEdge
                case 1
                    nodesPatch=[down_nodes1];
                case 2
                    nodesPatch =[right_nodes1];
                case 3
                    nodesPatch =[up_nodes1];
                case 4
                    nodesPatch =[left_nodes1];
            end
            
            for indexElem=1:length(elemsPatch)
                for indexNode=1:length(nodesPatch)
                    PHUTelem{indexPatch}(elemsPatch(indexElem)).solIndex( nodesPatch)=PHUTelem{neighborPatch}(elemsNeighbor(indexElem)).solIndex(nodesNeighbor);
                end
            end
            
            
        end
        
    end
    
    
    %disp('dealing with elem within patches')
    edge=patchInfo{indexPatch}(:,3);

    
    for i=1:length(edge)

        [~,elems] = sortEdgeNodesElem( PHUTelem{indexPatch}, edge(i), p, q);
        elems=cell2mat(elems);

        switch edge(i)
            case 1
                nodes=[down_nodes1,down_nodes2];
            case 2
                nodes =[right_nodes1,right_nodes2];
            case 3
                nodes =[up_nodes1,up_nodes2];
            case 4
                nodes =[left_nodes1,left_nodes2];
        end

        for indexElem=1:length(elems)

            %%%%%%%%assign solIndex From neighbor elem %%%%%%%%%
            neighbor=cell(1,4);
            neighbor{1}=PHUTelem{indexPatch}(elems(indexElem)).neighbor_down;
            neighbor{2}=PHUTelem{indexPatch}(elems(indexElem)).neighbor_up;
            neighbor{3}=PHUTelem{indexPatch}(elems(indexElem)).neighbor_left;
            neighbor{4}=PHUTelem{indexPatch}(elems(indexElem)).neighbor_right;
            
            for indexCell=1:4
                
                if ~isempty(neighbor{indexCell})
                    if neighbor{indexCell}<elems(indexElem)
                        
                        
                         switch indexCell
                            case 1 %  neighbor down
                                neighborNodes=up_nodes1;
                                elemNodes=down_nodes1;
                            case 2 % neighbor up
                                neighborNodes=down_nodes1;
                                elemNodes=up_nodes1;
                            case 3 % neighbor left
                                neighborNodes=right_nodes1;
                                elemNodes=left_nodes1;
                            case 4 % neighbor right
                                neighborNodes=left_nodes1;
                                elemNodes=right_nodes1;
                         end
                        
                        neighborElem=neighbor{indexCell};
                        for ii=1:length(neighborNodes)
                            if PHUTelem{indexPatch}(neighborElem).solIndex(neighborNodes(ii))~=0
                                PHUTelem{indexPatch}(elems(indexElem)).solIndex( elemNodes(ii))=PHUTelem{indexPatch}(neighborElem).solIndex(neighborNodes(ii));
                            end
                        end
                        
                    end
                    
                end
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for indexNode=1:length(nodes)
                if PHUTelem{indexPatch}(elems(indexElem)).solIndex(nodes(indexNode))==0
                    PHUTelem{indexPatch}(elems(indexElem)).solIndex(nodes(indexNode))=solIndexCount+1;
                    solIndexCount=solIndexCount+1;
                end
            end
            
        end
        
    end
    
    
    
end


end



