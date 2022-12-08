function [coefSol] =zipConformingC1_p2_BC_3PatchBiharmonic( PHUTelem,m,p,q )

right_nodes2 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
right_nodes1 = (p+1):(p+1):(p+1)*(q+1);

left_nodes1 = 1:(p+1):(1+(p+1)*q);
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

down_nodes1=1:(p+1);
down_nodes2=down_nodes1+(p+1);

up_nodes1=(p+1)*(q+1)-p:(p+1)*(q+1);
up_nodes2=up_nodes1-(p+1);
list=[];

indexPatch=1;
edge=2;
tempElem= sortEdgeElem( PHUTelem{indexPatch}, edge);
tempElem=cell2mat(tempElem);
for indexElem=tempElem
    if isempty(PHUTelem{indexPatch}(indexElem).neighbor_up)
        tempList=PHUTelem{indexPatch}(indexElem).solIndex([right_nodes1,right_nodes2]);
        list=[list,tempList];
    end
end

edge=4;
tempElem= sortEdgeElem( PHUTelem{indexPatch}, edge);
tempElem=cell2mat(tempElem);
for indexElem=tempElem
    if isempty(PHUTelem{indexPatch}(indexElem).neighbor_down)
        tempList=PHUTelem{indexPatch}(indexElem).solIndex([down_nodes1,down_nodes2]);
        list=[list,tempList];
    end
end


indexPatch=2;
edge=2;
tempElem= sortEdgeElem( PHUTelem{indexPatch}, edge);
tempElem=cell2mat(tempElem);
for indexElem=tempElem
    if isempty(PHUTelem{indexPatch}(indexElem).neighbor_down)
        tempList=PHUTelem{indexPatch}(indexElem).solIndex([right_nodes1,right_nodes2]);
        list=[list,tempList];
    end
end

edge=4;
tempElem= sortEdgeElem( PHUTelem{indexPatch}, edge);
tempElem=cell2mat(tempElem);
for indexElem=tempElem
    if isempty(PHUTelem{indexPatch}(indexElem).neighbor_up)
        tempList=PHUTelem{indexPatch}(indexElem).solIndex([up_nodes1,up_nodes2]);
        list=[list,tempList];
    end
end

indexPatch=3;
edge=2;
tempElem= sortEdgeElem( PHUTelem{indexPatch}, edge);
tempElem=cell2mat(tempElem);
for indexElem=tempElem
    if isempty(PHUTelem{indexPatch}(indexElem).neighbor_up)
        tempList=PHUTelem{indexPatch}(indexElem).solIndex([up_nodes1,up_nodes2]);
        list=[list,tempList];
    end
end

edge=4;
tempElem= sortEdgeElem( PHUTelem{indexPatch}, edge);
tempElem=cell2mat(tempElem);
for indexElem=tempElem
    if isempty(PHUTelem{indexPatch}(indexElem).neighbor_down)
        tempList=PHUTelem{indexPatch}(indexElem).solIndex([down_nodes1,down_nodes2]);
        list=[list,tempList];
    end
end

list=nonzeros(list);
list=unique(list);

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


%coefSol = nullMDS(m);
%numType2Basis=size(coefSol,2)

end

