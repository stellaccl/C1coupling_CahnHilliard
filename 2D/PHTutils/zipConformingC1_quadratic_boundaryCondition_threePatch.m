function [coefSol,numType2Basis] = zipConformingC1_quadratic_boundaryCondition_threePatch( PHUTelem,m,p,q )

right_nodes2 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
right_nodes1 = (p+1):(p+1):(p+1)*(q+1);

left_nodes1 = 1:(p+1):(1+(p+1)*q);
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

down_nodes1=1:(p+1);
down_nodes2=down_nodes1+(p+1);

up_nodes1=(p+1)*(q+1)-p:(p+1)*(q+1);
up_nodes2=up_nodes1-(p+1);


elemPatch3 = sortEdgeElem( PHUTelem{3}, 4);
elemPatch3=cell2mat(elemPatch3);
elem=[];

for elemIndex=1:length(elemPatch3)
    
    down=PHUTelem{3}(elemPatch3(elemIndex)).neighbor_down;
    
    if isempty(down)
        elem=[elemPatch3(elemIndex)];
    end
    
end

list=[];
for indexElem=1:length(elem)
    list=[list,PHUTelem{3}(elem(indexElem)).solIndex(left_nodes1)];
end
list=nonzeros(list);

m(:,list) = [];
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

numType2Basis=size(coefSol,2)

end

