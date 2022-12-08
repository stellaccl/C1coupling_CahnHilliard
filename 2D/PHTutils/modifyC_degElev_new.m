function [ PHUTelem ] = modifyC_degElev_new( PHUTelem,PHUTelemDE,coefSol,numType2Basis,patchEdgeInfo,boundaryInfo,type2Basis,p,q)
%compute modifiedC for degree elevation
%for 2 patches


right_nodes2 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
right_nodes1 = (p+1):(p+1):(p+1)*(q+1);

left_nodes1 = 1:(p+1):(1+(p+1)*q);
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

down_nodes1=1:(p+1);
down_nodes2=down_nodes1+(p+1);

up_nodes1=(p+1)*(q+1)-p:(p+1)*(q+1);
up_nodes2=up_nodes1-(p+1);

right_nodes3DE = (p+2)-2:(p+2):(p+2)*(q+2)-2;
right_nodes2DE = (p+2)-1:(p+2):(p+2)*(q+2)-1;
right_nodes1DE = (p+2):(p+2):(p+2)*(q+2);

left_nodes1DE = 1:(p+2):(1+(p+2)*(q+1));
left_nodes2DE = 2:(p+2):(1+(p+2)*(q+1))+1;
left_nodes3DE = 3:(p+2):(1+(p+2)*(q+1))+2;

down_nodes1DE=1:(p+2);
down_nodes2DE=down_nodes1DE+(p+2);
down_nodes3DE=down_nodes2DE+(p+2);

up_nodes1DE=(p+2)*(q+2)-(p+1):(p+2)*(q+2);
up_nodes2DE=up_nodes1DE-(p+2);
up_nodes3DE=up_nodes2DE-(p+2);

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
            
            if ismember(indexElem,boundaryElem)
                
                neighbor_up=PHUTelem{indexPatch}(indexElem).neighbor_up;
                neighbor_down=PHUTelem{indexPatch}(indexElem).neighbor_down;
                neighbor_left=PHUTelem{indexPatch}(indexElem).neighbor_left;
                neighbor_right=PHUTelem{indexPatch}(indexElem).neighbor_right;
                
                coefIndex=[];
                solIndex=[];
                
                for indexEdge=1:length(type2Edge)
                    
                    switch type2Edge(indexEdge)
                        
                        case 1
                            if isempty(neighbor_down)
                                %    disp('enter down')
                                coefIndex=[coefIndex,down_nodes2DE,down_nodes1DE];
                                solIndex=[solIndex,PHUTelemDE{indexPatch}(indexElem).solIndex(down_nodes2DE),PHUTelemDE{indexPatch}(indexElem).solIndex(down_nodes1DE)];
                                tempAdditionalBasisIndex=down_nodes3DE;
                            end
                            
                        case 2
                            if isempty(neighbor_right)
                                % disp('enter right')
                                coefIndex=[coefIndex,right_nodes2DE,right_nodes1DE];
                                solIndex=[solIndex,PHUTelemDE{indexPatch}(indexElem).solIndex(right_nodes2DE),PHUTelemDE{indexPatch}(indexElem).solIndex(right_nodes1DE)];
                                tempAdditionalBasisIndex=right_nodes3DE;
                            end
                            
                        case 3
                            if isempty(neighbor_up)
                                % disp('enter up')
                                coefIndex=[coefIndex,up_nodes2DE,up_nodes1DE];
                                solIndex=[solIndex,PHUTelemDE{indexPatch}(indexElem).solIndex(up_nodes2DE),PHUTelemDE{indexPatch}(indexElem).solIndex(up_nodes1DE)];
                                tempAdditionalBasisIndex=up_nodes3DE;
                            end
                            
                        case 4
                            if isempty(neighbor_left)
                                %   disp('enter left')
                                coefIndex=[coefIndex,left_nodes2DE,left_nodes1DE];
                                solIndex=[solIndex,PHUTelemDE{indexPatch}(indexElem).solIndex(left_nodes2DE),PHUTelemDE{indexPatch}(indexElem).solIndex(left_nodes1DE)];
                                tempAdditionalBasisIndex=left_nodes3DE;
                            end
                            
                    end
                    
                end
                
                [~,loc]=ismember(type2Basis,PHUTelem{indexPatch}(indexElem).nodesGlobal);
                numBasis=length(PHUTelem{indexPatch}(indexElem).nodesGlobal);
                
                type2NodesIndex=loc;
                
                tempC = PHUTelem{indexPatch}(indexElem).modifiedC;
                PHUTelem{indexPatch}(indexElem).modifiedC=zeros(numBasis,size(PHUTelemDE{indexPatch}(indexElem).C,2));
                PHUTelem{indexPatch}(indexElem).modifiedC(1:size(tempC,1),1:size(tempC,2)) = tempC;
                
                for ii=1:numType2Basis
                    localIndex =  type2NodesIndex(ii);
                    PHUTelem{indexPatch}(indexElem).modifiedC(localIndex,:) = 0;
                    PHUTelem{indexPatch}(indexElem).modifiedC(localIndex,coefIndex) = coefSol(solIndex,ii);
                end
                
                additionalBasisIndex=PHUTelem{indexPatch}(indexElem).additionalBasisIndex;
                PHUTelem{indexPatch}(indexElem).modifiedC(additionalBasisIndex,:) = PHUTelemDE{indexPatch}(indexElem).C(tempAdditionalBasisIndex,:);
                
            end
            
        end
    end
    
end
end

