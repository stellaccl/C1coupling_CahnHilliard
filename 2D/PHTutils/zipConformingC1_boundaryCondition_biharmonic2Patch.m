function [ coefSol ] = zipConformingC1_boundaryCondition_biharmonic2Patch( PHUTelem,m,p,q )
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
elemPatch1 = sortEdgeElem( PHUTelem{1}, 2);
elemPatch2 = sortEdgeElem( PHUTelem{2}, 4);
elemPatch1=cell2mat(elemPatch1);
elemPatch2=cell2mat(elemPatch2);

elemUpPatch1=elemPatch1(end);
elemDownPatch1=elemPatch1(1);

elemUpPatch2=elemPatch2(end);
elemDownPatch2=elemPatch2(1);

list=[];
% list=[list,PHUTelem{1}(elemUpPatch1).solIndex(up_nodes1)];
% list=[list,PHUTelem{2}(elemUpPatch2).solIndex(up_nodes1)];
% list=[list,PHUTelem{1}(elemDownPatch1).solIndex(down_nodes1)];
% list=[list,PHUTelem{2}(elemDownPatch2).solIndex(down_nodes1)];

list=[list,PHUTelem{1}(elemUpPatch1).solIndex(up_nodes1),PHUTelem{1}(elemUpPatch1).solIndex(up_nodes2)];
list=[list,PHUTelem{2}(elemUpPatch2).solIndex(up_nodes1),PHUTelem{2}(elemUpPatch2).solIndex(up_nodes2)];
list=[list,PHUTelem{1}(elemDownPatch1).solIndex(down_nodes1),PHUTelem{1}(elemDownPatch1).solIndex(down_nodes2)];
list=[list,PHUTelem{2}(elemDownPatch2).solIndex(down_nodes1),PHUTelem{2}(elemDownPatch2).solIndex(down_nodes2)];

list=nonzeros(list);
list=unique(list);

m(:,list) = [];
% tempCoefSol = nullMDS(m);
tempCoefSol = null(m);
%add rows of zeros
coefSol=zeros(size(tempCoefSol,1)+length(list),size(tempCoefSol,2));
count=0;
for indexRow=1:size(coefSol,1)
    if  ~ismember(indexRow,list)
        count=count+1;
        coefSol(indexRow,:)=tempCoefSol(count,:);
    end
end

