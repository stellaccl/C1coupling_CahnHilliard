function [ coefSol ] = zipConformingC1_boundaryCondition_bilinearCahnHilliard( PHUTelem,m,p,q )
%apply boundary condition to m

right_nodes2 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
right_nodes1 = (p+1):(p+1):(p+1)*(q+1);

left_nodes1 = 1:(p+1):(1+(p+1)*q);
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

down_nodes1=1:(p+1);
down_nodes2=down_nodes1+(p+1);

up_nodes1=(p+1)*(q+1)-p:(p+1)*(q+1);
up_nodes2=up_nodes1-(p+1);

tempM=zeros(size(m,1)+3,size(m,2)); %add 3 constraints
tempM(1:size(m,1),1:size(m,2))=m;

elemA=sortEdgeElem(PHUTelem{1},2);
elemB=sortEdgeElem(PHUTelem{2},4);

% elemA_up=elemA(end);
% elemA_down=elemA(1);
% 
% elemB_up=elemB(end);
% elemB_down=elemB(1);

elemA=cell2mat(elemA);
elemB=cell2mat(elemB);

list=zeros(3,3);
list(1,:)=[PHUTelem{1}(elemA(1)).solIndex([3,7]),PHUTelem{1}(elemA(end)).solIndex(11)];
list(2,:)=[PHUTelem{1}(elemA(1)).solIndex([4,8]),PHUTelem{1}(elemA(end)).solIndex(12)];
list(3,:)=[PHUTelem{2}(elemB(1)).solIndex([2,6]),PHUTelem{2}(elemB(end)).solIndex(10)];

for i=0:2
    tempM(end-i,list(i+1,:))=[2,-1,-1];
end

coefSol = nullMDS(tempM);

%======================= apply boundary condition =========================
% boundaryPatch=11:18;
% list=[];
% for indexPatch=1:length(boundaryPatch)
%     
%     %select element at the left side
%     elem=sortEdgeElem(PHUTelem{boundaryPatch(indexPatch)},4);
%     elem=cell2mat(elem);
%     
%     for indexElem=1:length(elem)
%         list=[list,PHUTelem{boundaryPatch(indexPatch)}(elem(indexElem)).solIndex(left_nodes1)];
%     end
%     
% end
% 
% list=nonzeros(list);
% list=unique(list);
% 
% m(:,list) = [];
% svd(m)
% tempCoefSol = null(m);
% 
% %tempCoefSol = nullMDS(m);
% %add rows of zeros
% coefSol=zeros(size(tempCoefSol,1)+length(list),size(tempCoefSol,2));
% count=0;
% for indexRow=1:size(coefSol,1)
%     if  ~ismember(indexRow,list)
%         count=count+1;
%         coefSol(indexRow,:)=tempCoefSol(count,:);
%     end
% end


