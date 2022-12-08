function [ PHUTelem ] = modifyC_3( PHUTelem,coefSol,numType2Basis,patchEdgeInfo,boundaryInfo,type2Basis,p,q)
%compute modifiedC for five patch example

right_nodes2 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
right_nodes1 = (p+1):(p+1):(p+1)*(q+1);

left_nodes1 = 1:(p+1):(1+(p+1)*q);
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

down_nodes1=1:(p+1);
down_nodes2=down_nodes1+(p+1);

up_nodes1=(p+1)*(q+1)-p:(p+1)*(q+1);
up_nodes2=up_nodes1-(p+1);


for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        PHUTelem{indexPatch}(indexElem).modifiedC=PHUTelem{indexPatch}(indexElem).C;
    end
end

for indexPatch=1:length(PHUTelem)
    
    type2Edge=patchEdgeInfo{indexPatch};
    
    boundaryElem=boundaryInfo{indexPatch}(:,5);
    boundaryElem=unique(boundaryElem);
    
    for indexElem=1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(indexElem).children)
            coefIndex=[];
            solIndex=[];
            if ismember(indexElem,boundaryElem)
                
                neighbor_up=PHUTelem{indexPatch}(indexElem).neighbor_up;
                neighbor_down=PHUTelem{indexPatch}(indexElem).neighbor_down;
                neighbor_left=PHUTelem{indexPatch}(indexElem).neighbor_left;
                neighbor_right=PHUTelem{indexPatch}(indexElem).neighbor_right;
                
                for indexEdge=1:length(type2Edge)
                    
                    switch type2Edge(indexEdge)
                        
                        case 1
                            if isempty(neighbor_down)
                                %    disp('enter down')
                                coefIndex=[coefIndex,down_nodes2,down_nodes1];
                                solIndex=[solIndex,PHUTelem{indexPatch}(indexElem).solIndex(down_nodes2),PHUTelem{indexPatch}(indexElem).solIndex(down_nodes1)];
                            end
                            
                        case 2
                            if isempty(neighbor_right)
                                % disp('enter right')
                                coefIndex=[coefIndex,right_nodes2,right_nodes1];
                                solIndex=[solIndex,PHUTelem{indexPatch}(indexElem).solIndex(right_nodes2),PHUTelem{indexPatch}(indexElem).solIndex(right_nodes1)];
                            end
                            
                        case 3
                            if isempty(neighbor_up)
                                % disp('enter up')
                                coefIndex=[coefIndex,up_nodes2,up_nodes1];
                                solIndex=[solIndex,PHUTelem{indexPatch}(indexElem).solIndex(up_nodes2),PHUTelem{indexPatch}(indexElem).solIndex(up_nodes1)];
                            end
                            
                        case 4
                            if isempty(neighbor_left)
                                %   disp('enter left')
                                coefIndex=[coefIndex,left_nodes2,left_nodes1];
                                solIndex=[solIndex,PHUTelem{indexPatch}(indexElem).solIndex(left_nodes2),PHUTelem{indexPatch}(indexElem).solIndex(left_nodes1)];
                            end
                            
                    end
                    
                end
                
                [~,loc]=ismember(type2Basis,PHUTelem{indexPatch}(indexElem).nodesGlobal);
                type2NodesIndex=loc;
                
                for ii=1:numType2Basis
                    localIndex =  type2NodesIndex(ii);
                    PHUTelem{indexPatch}(indexElem).modifiedC(localIndex,:) = 0;
                    PHUTelem{indexPatch}(indexElem).modifiedC(localIndex,coefIndex) = coefSol(solIndex,ii);
                end
                
            end
        end
        
    end
    
end


end

