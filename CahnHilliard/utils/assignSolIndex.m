function [ PHUTelem,numConstraints ] = assignSolIndex( PHUTelem,patchBoundaries,p,q )
numBoundaries = size(patchBoundaries,1);
patchesSeen = [];
for boundaryIndex = 1:numBoundaries
    
    patchAList = patchBoundaries{boundaryIndex,1};
    patchB = patchBoundaries{boundaryIndex,2};
    edgeAList = patchBoundaries{boundaryIndex,3};
    edgeBList = patchBoundaries{boundaryIndex,4};
    
    for indexPatch=1:length(patchAList)
        
        patchA = patchAList(indexPatch);
        edgeA = edgeAList(indexPatch);
        edgeB = edgeBList(indexPatch);
        patchB = patchB(indexPatch);
        
        ElemA = sortEdgeElem( PHUTelem{patchA}, edgeA);
        ElemB = sortEdgeElem( PHUTelem{patchB}, edgeB);
    end
end

right_nodes1 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
right_nodes2 = (p+1):(p+1):(p+1)*(q+1);

left_nodes1 = 1:(p+1):(1+(p+1)*q);
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

down_nodes1= 1:p+1;
down_nodes2=down_nodes1+(p+1);

up_nodes2=(p+1)*(q+1)-p:(p+1)*(q+1);
up_nodes1=up_nodes2-(p+1);
count=0;
for i=1:length(ElemA{1})
    
    PHUTelem{1}(ElemA{1}(i)).solIndex=zeros(1,16);
    PHUTelem{2}(ElemB{1}(i)).solIndex=zeros(1,16);
    for j=1:p+1
        count=count+1;
        PHUTelem{1}(ElemA{1}(i)).solIndex(right_nodes1(j))=count;
        count=count+1;
        PHUTelem{1}(ElemA{1}(i)).solIndex(right_nodes2(j))=count;
        PHUTelem{2}(ElemB{1}(i)).solIndex(left_nodes1(j))=count;
        count=count+1;
        PHUTelem{2}(ElemB{1}(i)).solIndex(left_nodes2(j))=count;   
    end
   
    count=count-3;
    
%     if i==length(ElemA{1})
%         PHUTelem{1}(ElemA{1}(i)).solIndex([11,12,15,16])=PHUTelem{1}(ElemA{1}(1)).solIndex([7,8,3,4]);
%         PHUTelem{2}(ElemB{1}(i)).solIndex([9,10,13,14])=PHUTelem{2}(ElemB{1}(1)).solIndex([5,6,1,2]);
%         numConstraints=PHUTelem{2}(ElemB{1}(i)).solIndex(6); %numConstrain in matrix m
%     end

    if i==length(ElemA{1})
        PHUTelem{1}(ElemA{1}(i)).solIndex([15,16])=PHUTelem{1}(ElemA{1}(1)).solIndex([3,4]);
        PHUTelem{2}(ElemB{1}(i)).solIndex([13,14])=PHUTelem{2}(ElemB{1}(1)).solIndex([1,2]);
        numConstraints=PHUTelem{2}(ElemB{1}(i)).solIndex(10); %numConstrain in matrix m
    end


end

end

