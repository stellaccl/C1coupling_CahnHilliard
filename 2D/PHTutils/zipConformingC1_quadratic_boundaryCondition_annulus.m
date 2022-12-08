function [coefSol] = zipConformingC1_quadratic_boundaryCondition_annulus( PHUTelem,m,p,q )

right_nodes2 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
right_nodes1 = (p+1):(p+1):(p+1)*(q+1);

left_nodes1 = 1:(p+1):(1+(p+1)*q);
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

down_nodes1=1:(p+1);
down_nodes2=down_nodes1+(p+1);

up_nodes1=(p+1)*(q+1)-p:(p+1)*(q+1);
up_nodes2=up_nodes1-(p+1);

list=[];
for indexPatch=1:length(PHUTelem)
    
    %select element at the left side
    elem_left=sortEdgeElem(PHUTelem{indexPatch},4);
    elem_right=sortEdgeElem(PHUTelem{indexPatch},2);
    
    elem_left=cell2mat(elem_left);
    elem_right=cell2mat(elem_right);
    
    for indexElem=1:length(elem_left)
        list=[list,PHUTelem{indexPatch}(elem_left(indexElem)).solIndex(left_nodes1),PHUTelem{indexPatch}(elem_left(indexElem)).solIndex(left_nodes2)];
    end
    
    for indexElem=1:length(elem_right)
        list=[list,PHUTelem{indexPatch}(elem_right(indexElem)).solIndex(right_nodes1),PHUTelem{indexPatch}(elem_right(indexElem)).solIndex(right_nodes2)];
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

