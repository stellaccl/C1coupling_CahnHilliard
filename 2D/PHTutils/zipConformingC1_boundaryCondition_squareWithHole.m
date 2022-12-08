function [coefSol] = zipConformingC1_boundaryCondition_squareWithHole( PHUTelem,m,p,q )

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

list=[];
for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
      if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_right)
          list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesRight)];
      end
    end
end

for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
      if  isempty(PHUTelem{indexPatch}(indexElem).neighbor_left)
          list=[list,PHUTelem{indexPatch}(indexElem).solIndex(nodesLeft)];
      end
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

