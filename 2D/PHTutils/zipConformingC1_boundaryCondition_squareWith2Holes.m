function [coefSol] = zipConformingC1_boundaryCondition_squareWith2Holes( PHUTelem,m,p,q )

right_nodes2 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
right_nodes1 = (p+1):(p+1):(p+1)*(q+1);

left_nodes1 = 1:(p+1):(1+(p+1)*q);
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

down_nodes1=1:(p+1);
down_nodes2=down_nodes1+(p+1);

up_nodes1=(p+1)*(q+1)-p:(p+1)*(q+1);
up_nodes2=up_nodes1-(p+1);

nodesRight=[right_nodes1, right_nodes2];
nodesLeft=[left_nodes1, left_nodes2];
nodesUp=[up_nodes1, up_nodes2];
nodesDown=[down_nodes1, down_nodes2];

list=[];

indexPatch=1;
for indexElem=1:length(PHUTelem{indexPatch})
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_right)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesRight)];
    end
    
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_left)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesLeft)];
    end
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indexPatch=2;
for indexElem=1:length(PHUTelem{indexPatch})

    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_down)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesDown)];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indexPatch=3;
for indexElem=1:length(PHUTelem{indexPatch})

    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_up)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesUp)];
    end
    
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_down)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesDown)];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indexPatch=4;
for indexElem=1:length(PHUTelem{indexPatch})

    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_left)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesLeft)];
    end

    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_down)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesDown)];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indexPatch=5;
for indexElem=1:length(PHUTelem{indexPatch})
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_right)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesRight)];
    end
    
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_left)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesLeft)];
    end
   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indexPatch=6;
for indexElem=1:length(PHUTelem{indexPatch})
  
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_left)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesLeft)];
    end
    
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_up)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesUp)];
    end
 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indexPatch=7;
for indexElem=1:length(PHUTelem{indexPatch})
   
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_up)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesUp)];
    end
    
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_down)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesDown)];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indexPatch=8;
for indexElem=1:length(PHUTelem{indexPatch})

    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_up)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesUp)];
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indexPatch=9;
for indexElem=1:length(PHUTelem{indexPatch})

    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_up)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesUp)];
    end
    
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_down)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesDown)];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indexPatch=10;
for indexElem=1:length(PHUTelem{indexPatch})
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_right)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesRight)];
    end

    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_up)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesUp)];
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indexPatch=11;
for indexElem=1:length(PHUTelem{indexPatch})
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_right)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesRight)];
    end
    
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_left)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesLeft)];
    end
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indexPatch=12;
for indexElem=1:length(PHUTelem{indexPatch})
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_right)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesRight)];
    end

    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_down)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesDown)];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indexPatch=13;
for indexElem=1:length(PHUTelem{indexPatch})
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_up)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesUp)];
    end
    
    if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_down)
        list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesDown)];
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

