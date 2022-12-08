function [PHUTelem,m] =zipConformingC1_shell_method2_multiElemGIFTmesh(PHUTelem,GIFTmesh,patchBoundaries,solIndexCount,p,q)

[boundaryInfo,segInfo,~,~,patchNeighborInfo,patchEdgeInfo] = createBoundaryInfo3(PHUTelem, patchBoundaries );

numSeg=size(segInfo,1);

numBoundaries = size(patchBoundaries,1);

numColumns=solIndexCount;
numPts=ceil(numColumns/numSeg);
if numPts<p+1
    numPts=p+1;
end

m=zeros(numPts*numSeg*3,numColumns);
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
        
        [m,i] =computeMatrixM_shell_method2_multiElemGIFTmesh(PHUTelem,GIFTmesh,patchInfo,m,i,p,q,numPts);
        
    end
    
end


end

