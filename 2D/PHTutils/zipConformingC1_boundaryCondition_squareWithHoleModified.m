function [coefSol] = zipConformingC1_boundaryCondition_squareWithHoleModified( PHUTelem,m,p,q )

right_nodes2 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
right_nodes1 = (p+1):(p+1):(p+1)*(q+1);

left_nodes1 = 1:(p+1):(1+(p+1)*q);
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

down_nodes1=1:(p+1);
down_nodes2=down_nodes1+(p+1);

up_nodes1=(p+1)*(q+1)-p:(p+1)*(q+1);
up_nodes2=up_nodes1-(p+1);

leftNodes=[left_nodes1, left_nodes2];
rightNodes=[right_nodes1, right_nodes2];
upNodes=[up_nodes1, up_nodes2];
downNodes=[down_nodes1, down_nodes2];

list=[];
for indexPatch=1
    for indexElem=1:length(PHUTelem{indexPatch})
        if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_up)
            list=[list,PHUTelem{indexPatch}(indexElem).solIndex(upNodes)];
        end
        
        if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_down)
            list=[list,PHUTelem{indexPatch}(indexElem).solIndex(downNodes)];
        end
        
    end
end

indexPatch=2;
for indexElem=1:length(PHUTelem{indexPatch})
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_up)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(upNodes)];
    end
    
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_left)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(leftNodes)];
    end
end

indexPatch=3;
for indexElem=1:length(PHUTelem{indexPatch})
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_right)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(rightNodes)];
    end
    
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_left)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(leftNodes)];
    end
end


indexPatch=4;
for indexElem=1:length(PHUTelem{indexPatch})
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_down)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(downNodes)];
    end
    
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_left)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(leftNodes)];
    end
end

indexPatch=5;
for indexElem=1:length(PHUTelem{indexPatch})
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_down)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(downNodes)];
    end
    
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_up)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(upNodes)];
    end
end

indexPatch=6;
for indexElem=1:length(PHUTelem{indexPatch})
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_down)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(downNodes)];
    end
    
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_right)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(rightNodes)];
    end
end


indexPatch=7;
for indexElem=1:length(PHUTelem{indexPatch})
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_left)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(leftNodes)];
    end
    
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_right)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(rightNodes)];
    end
end

indexPatch=8;
for indexElem=1:length(PHUTelem{indexPatch})
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_up)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(upNodes)];
    end
    
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_right)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(rightNodes)];
    end
end

list=nonzeros(list);
list=unique(list);

m(:,list) = [];

%tempCoefSol = null(m);

tempCoefSol = nullMDS(m);
%add rows of zeros
coefSol=zeros(size(tempCoefSol,1)+length(list),size(tempCoefSol,2));
count=0;
for indexRow=1:size(coefSol,1)
    if  ~ismember(indexRow,list)
        count=count+1;
        coefSol(indexRow,:)=tempCoefSol(count,:);
    end
end

end

