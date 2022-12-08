function [PHUTelem,sizeBasis,numType2Basis,type2Basis] =zipConformingC1_degElev(PHUTelem,PHUTelemDE,GIFTmesh,sizeBasis,patchBoundaries,solIndexCount,p,q)
%zipConforming without boundary condition
disp('in zipConforming C1_degElev test')
[boundaryInfo,segInfo,subSegInfo,numSeg,patchNeighborInfo,patchEdgeInfo] = createBoundaryInfo3(PHUTelem, patchBoundaries );

numBoundaries = size(patchBoundaries,1);

numColumns=solIndexCount;
numPts=ceil(numColumns/numSeg);
%numPts correspond to degree elevated elements
if numPts<p+2
    numPts=p+2;
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
        
        [m,i] = computeMatrixM(PHUTelemDE,GIFTmesh,patchInfo,m,i,p+1,q+1,numPts);
        
    end
    
end

coefSol = nullMDS(m);
numType2Basis=size(coefSol,2)
type2basisNodes=sizeBasis+1:sizeBasis+numType2Basis;
[PHUTelem,type2Basis,sizeBasis] = assignNewNodesGlobal_degElev(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,type2basisNodes,p,q);

[PHUTelem] = modifyC_degElev_new(PHUTelem,PHUTelemDE,coefSol,numType2Basis,patchEdgeInfo,boundaryInfo,type2Basis,p,q);

end

