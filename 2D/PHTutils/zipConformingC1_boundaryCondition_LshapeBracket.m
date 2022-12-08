function [ coefSol ] = zipConformingC1_boundaryCondition_LshapeBracket( PHUTelem,m,p,q )
%apply boundary condition to m

right_nodes2 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
right_nodes1 = (p+1):(p+1):(p+1)*(q+1);

left_nodes1 = 1:(p+1):(1+(p+1)*q);
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

down_nodes1=1:(p+1);
down_nodes2=down_nodes1+(p+1);

up_nodes1=(p+1)*(q+1)-p:(p+1)*(q+1);
up_nodes2=up_nodes1-(p+1);

%======================= apply boundary condition =========================
boundaryPatch=11:18;
list=[];
for indexPatch=1:length(boundaryPatch)
    
    %select element at the left side
    elem=sortEdgeElem(PHUTelem{boundaryPatch(indexPatch)},4);
    elem=cell2mat(elem);
    
    for indexElem=1:length(elem)
        list=[list,PHUTelem{boundaryPatch(indexPatch)}(elem(indexElem)).solIndex(left_nodes1)];
    end
    
end

list=nonzeros(list);
list=unique(list);

m(:,list) = [];
svd(m)
tempCoefSol = null(m);

%tempCoefSol = nullMDS(m);
%add rows of zeros
coefSol=zeros(size(tempCoefSol,1)+length(list),size(tempCoefSol,2));
count=0;
for indexRow=1:size(coefSol,1)
    if  ~ismember(indexRow,list)
        count=count+1;
        coefSol(indexRow,:)=tempCoefSol(count,:);
    end
end


