function [PHUTelem,dimBasis,numType2Basis] =zipConforming_c1_arbitraryDegreePB3(PHUTelem,GIFTmesh,dimBasis, patchBoundaries, p, q,numConstraints)
%assign new global nodes
%use assigned sol index in PHUTelem
numPatches = length(PHUTelem);
numBoundaries = size(patchBoundaries,1);

right_nodes1 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
right_nodes2 = (p+1):(p+1):(p+1)*(q+1);

left_nodes1 = 1:(p+1):(1+(p+1)*q);
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

for boundaryIndex = 1:numBoundaries
    patchA = patchBoundaries(boundaryIndex,1);
    patchB = patchBoundaries(boundaryIndex,2);
    edgeA = patchBoundaries(boundaryIndex,3);
    edgeB = patchBoundaries(boundaryIndex,4);
    
    edgeA=cell2mat(edgeA);
    edgeB=cell2mat(edgeB);
    
    patchA=cell2mat(patchA);
    patchB=cell2mat(patchB);
    
    [elementsA] = sortEdgeElem( PHUTelem{patchA}, edgeA);
    [elementsB] = sortEdgeElem( PHUTelem{patchB}, edgeB);
    
    elementsA=elementsA{1};
    elementsB=elementsB{1};
    
    
    % [m] = computeMatrixM2_arbitraryDegree( PHUTelem,GIFTmesh,elementsA, elementsB,p,q );
    [m] = computeMatrixM2_arbitraryDegreePB2( PHUTelem,GIFTmesh,elementsA, elementsB,p,q,numConstraints );

    list=zeros(3,3);
    list(1,:)=[PHUTelem{1}(elementsA(1)).solIndex([3,7]),PHUTelem{1}(elementsA(end)).solIndex(11)];
    list(2,:)=[PHUTelem{1}(elementsA(1)).solIndex([4,8]),PHUTelem{1}(elementsA(end)).solIndex(12)];
    list(3,:)=[PHUTelem{2}(elementsB(1)).solIndex([2,6]),PHUTelem{2}(elementsB(end)).solIndex(10)];
    
    for i=0:2
        m(end-i,list(i+1,:))=[2,-1,-1];
    end

%     list1=[1,4,70];
%     list2=[2,5,71];
%     list3=[3,6,72];
%     
%     m(end-2,list1)=[2,-1,-1];
%     m(end-1,list2)=[2,-1,-1];
%     m(end,list3)=[2,-1,-1];
    
    %coefSol=null(m);
    coefSol = nullMDS(m);
    
    numType2Basis=size(coefSol,2);
    
    type2basisNodes=dimBasis+1:dimBasis+numType2Basis;
    removedNodes=[];
    
    %modify the nodesGlobal indices in patchA
    for indexElem = elementsA
        modifyBasisIndex=[right_nodes1,right_nodes2(1:end-1),right_nodes2(end):right_nodes2(end)+numType2Basis-(length(right_nodes2)*2)];
        
        removedNodes=[removedNodes,PHUTelem{patchA}(indexElem).nodesGlobal(modifyBasisIndex(1:length(right_nodes2)*2))];
        PHUTelem{patchA}(indexElem).nodesGlobal(modifyBasisIndex)=type2basisNodes;
        PHUTelem{patchA}(indexElem).type2nodes=modifyBasisIndex;
    end
    
    %modify the nodesGlobal indices in patchB
    for indexElem = elementsB
        modifyBasisIndex=[left_nodes1,left_nodes2,right_nodes2(end)+1:right_nodes2(end)+1+numType2Basis-(length(right_nodes2)*2+1)];
        removedNodes=[removedNodes,PHUTelem{patchB}(indexElem).nodesGlobal(modifyBasisIndex(1:length(right_nodes2)*2))];
        PHUTelem{patchB}(indexElem).nodesGlobal(modifyBasisIndex)=type2basisNodes;
        PHUTelem{patchB}(indexElem).type2nodes=modifyBasisIndex;
    end
    
    removedNodes=unique(removedNodes);
    
end

%re-assign global nodes
temp_nodes=1:dimBasis+numType2Basis;
replacePattern = temp_nodes;
temp_nodes(removedNodes)=[];
dimBasis=length(temp_nodes);
replacePattern(temp_nodes) = 1:dimBasis;

type2basis=replacePattern(type2basisNodes)

for indexPatch=1:numPatches
    for indexElem=1:length(PHUTelem{indexPatch})
        if isempty(PHUTelem{indexPatch}(indexElem).children)
            PHUTelem{indexPatch}(indexElem).nodesGlobal = replacePattern(PHUTelem{indexPatch}(indexElem).nodesGlobal);
        end
    end
end

%modify the Bezier coeficients to account for the C1 basis
%[PHUTelem] = modifyC_arbitraryDegree( PHUTelem,coefSol,numType2Basis,patchBoundaries,p,q);
%[PHUTelem] = modifyC_arbitraryDegreePB2( PHUTelem,coefSol,numType2Basis,patchBoundaries,p,q);
[PHUTelem] = modifyC_arbitraryDegreePB3( PHUTelem,coefSol,numType2Basis,patchBoundaries,p,q);




