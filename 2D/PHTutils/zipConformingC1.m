function [PHUTelem,sizeBasis,numType2Basis,type2Basis] =zipConformingC1(PHUTelem,GIFTmesh,sizeBasis,patchBoundaries,solIndexCount,p,q)
%zipConforming without boundary condition

right_nodes2 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
right_nodes1 = (p+1):(p+1):(p+1)*(q+1);

left_nodes1 = 1:(p+1):(1+(p+1)*q);
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

down_nodes1=1:(p+1);
down_nodes2=down_nodes1+(p+1);

up_nodes1=(p+1)*(q+1)-p:(p+1)*(q+1);
up_nodes2=up_nodes1-(p+1);

[boundaryInfo,segInfo,subSegInfo,numSeg,patchNeighborInfo,patchEdgeInfo] = createBoundaryInfo3(PHUTelem, patchBoundaries );

numBoundaries = size(patchBoundaries,1);

numColumns=solIndexCount;
numPts=ceil(numColumns/numSeg);
if numPts<p+1
    numPts=p+1;
end

m=zeros(numPts*numSeg,numColumns);
i=0;

for boundaryIndex = 1:numBoundaries
    patchAList = patchBoundaries(boundaryIndex,1);
    patchB = patchBoundaries(boundaryIndex,2);
    edgeAList = patchBoundaries(boundaryIndex,3);
    edgeBList = patchBoundaries(boundaryIndex,4);
    
    edgeAList=cell2mat(edgeAList);
    edgeBList=cell2mat(edgeBList);
    
    patchAList=cell2mat(patchAList);
    patchB=cell2mat(patchB);
    
    for indexPatch=1:length(patchAList)
        patchA = patchAList(indexPatch);
        edgeA = edgeAList(indexPatch);
        edgeB = edgeBList(indexPatch);
        
        [elemA] = sortEdgeElem( PHUTelem{patchA}, edgeA);
        [elemB] = sortEdgeElem( PHUTelem{patchB}, edgeB);
        
        elemA=elemA{1};
        elemB=elemB{1};
        
        patchInfo.patchA=patchA;
        patchInfo.patchB=patchB;
        patchInfo.edgeA=edgeA;
        patchInfo.edgeB=edgeB;
        patchInfo.elemA=elemA;
        patchInfo.elemB=elemB;

        [m,i] = computeMatrixM(PHUTelem,GIFTmesh,patchInfo,m,i,p,q,numPts);
    
    end
    
end

%%%%%%%%% boundary condition for circle biharmonic
% list=[];
% for indexPatch=1:4
%     
%     elem = sortEdgeElem( PHUTelem{indexPatch}, 1);
%     elem=cell2mat(elem);
%     
%     for indexElem=elem
%         list=[list,PHUTelem{indexPatch}(indexElem).solIndex(down_nodes1),PHUTelem{indexPatch}(indexElem).solIndex(down_nodes2)];
%     end
%     
% end
% list=nonzeros(list);
% list=unique(list);
% 
% 
% m(:,list) = [];
% tempCoefSol = nullMDS(m);
% 
% %add rows of zeros
% coefSol=zeros(size(tempCoefSol,1)+length(list),size(tempCoefSol,2));
% count=0;
% for indexRow=1:size(coefSol,1)
%     if  ~ismember(indexRow,list)
%         count=count+1;
%         coefSol(indexRow,:)=tempCoefSol(count,:);
%     end
% end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


coefSol = nullMDS(m);

numType2Basis=size(coefSol,2);


type2basisNodes=sizeBasis+1:sizeBasis+numType2Basis;
[PHUTelem ,type2Basis,sizeBasis] = assignNewNodesGlobal(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,type2basisNodes,p,q);
[PHUTelem] = modifyC(PHUTelem,coefSol,numType2Basis,patchEdgeInfo,boundaryInfo,type2Basis,p,q);

end

