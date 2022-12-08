function [ PHUTelem] = reassignNodesPR2( PHUTelem,numElem,p)
%change basis nodes

row=zeros((numElem-1)*2+(p+1),numElem*(p+1));
initial=1:p+1;
list=initial;
for i=1:numElem
    row(1,list)=initial+(i-1)*2;
    list=list+(p+1);
end
row(1,end-1)=row(1,1);
row(1,end)=row(1,2);

for i=2:size(row,1)-2
    row(i,:)=row(i-1,:)+size(row,1)-2;
end
row(end-1,:)=row(1,:);
row(end,:)=row(2,:);

count=0;

for rowIndex=1:numElem:numElem^2-(numElem-1)
    
    if rowIndex==1
        count=count+1;
    else
        count=count-1;
    end
    
    indexList=1:p+1;
    for i=0:numElem-1
        PHUTelem(rowIndex+i).nodes(1:4)=row(count,indexList);
        indexList=indexList+(p+1);
    end
    
    count=count+1;
    indexList=1:p+1;
    for i=0:numElem-1
        PHUTelem(rowIndex+i).nodes(5:8)=row(count,indexList);
        indexList=indexList+(p+1);
    end
    
    count=count+1;
    indexList=1:p+1;
    for i=0:numElem-1
        PHUTelem(rowIndex+i).nodes(9:12)=row(count,indexList);
        indexList=indexList+(p+1);
    end
    
    count=count+1;
    indexList=1:p+1;
    for i=0:numElem-1
        PHUTelem(rowIndex+i).nodes(13:16)=row(count,indexList);
        indexList=indexList+(p+1);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for indexElem=1:length(PHUTelem)
    
    PHUTelem(indexElem).nodesGlobal=PHUTelem(indexElem).nodes;
    
end


end

