function  [ PHUTelem,solIndexCount] = assignSolIndex2D_quadratic( PHUTelem,patchBoundaries,sizeBasis,p,q)
%assign solIndex
%support cross insertion
%assign remainNodes , the basis function not involving in C1 couplign and remaining in the solution space
% remainNodes : use for impose boundary condition in imposeDirichlet 

numBoundaries = size(patchBoundaries,1);
numPatches = length(PHUTelem);

%create/set nodesGlobal entries in all patches to be equal to local nodes
%entries
for patchIndex = 1:numPatches
    for elemIndex = 1:length(PHUTelem{patchIndex})
        PHUTelem{patchIndex}(elemIndex).solIndex = zeros(1,(p+1)*(q+1));
        PHUTelem{patchIndex}(elemIndex).remainNodes = PHUTelem{patchIndex}(elemIndex).nodesGlobal;
    end
end

right_nodes2 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
right_nodes1 = (p+1):(p+1):(p+1)*(q+1);

left_nodes1 = 1:(p+1):(1+(p+1)*q);
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

down_nodes1=1:(p+1);
down_nodes2=down_nodes1+(p+1);

up_nodes1=(p+1)*(q+1)-p:(p+1)*(q+1);
up_nodes2=up_nodes1-(p+1);

solIndexNodes=[];
for boundaryIndex = 1:numBoundaries
    
    patchAList = patchBoundaries{boundaryIndex,1};
    patchB = patchBoundaries{boundaryIndex,2};
    edgeAList = patchBoundaries{boundaryIndex,3};
    edgeBList = patchBoundaries{boundaryIndex,4};
    
    for indexPatch=1:length(patchAList)
        
        patchA = patchAList(indexPatch);
        edgeA = edgeAList(indexPatch);
        edgeB = edgeBList(indexPatch);
        
        elemA = sortEdgeElem( PHUTelem{patchA}, edgeA);
        elemB = sortEdgeElem( PHUTelem{patchB}, edgeB);
        
        elemA=cell2mat(elemA);
        elemB=cell2mat(elemB);
        
        switch edgeA
            case 1
                type2NodesA=[down_nodes1,down_nodes2];
            case 2
                type2NodesA=[right_nodes1,right_nodes2];
            case 3
                type2NodesA=[up_nodes1,up_nodes2];
            case 4
                type2NodesA=[left_nodes1,left_nodes2];
        end
        
        switch edgeB
            case 1
                type2NodesB=[down_nodes1,down_nodes2];
            case 2
                type2NodesB=[right_nodes1,right_nodes2];
            case 3
                type2NodesB=[up_nodes1,up_nodes2];
            case 4
                type2NodesB=[left_nodes1,left_nodes2];
        end
        
        for indexElem=1:length(elemA)
            solIndexNodes=[solIndexNodes,PHUTelem{patchA}(elemA(indexElem)).nodesGlobal(type2NodesA)];
            solIndexNodes=[solIndexNodes, PHUTelem{patchB}(elemB(indexElem)).nodesGlobal(type2NodesB)];
            
            PHUTelem{patchA}(elemA(indexElem)).remainNodes(type2NodesA)=0;
            PHUTelem{patchB}(elemB(indexElem)).remainNodes(type2NodesB)=0;
            
        end
        
    end
    
end
solIndexNodes=unique(solIndexNodes);
replacePattern=1:sizeBasis;
zeroSolIndex=setdiff(replacePattern,solIndexNodes);
replacePattern(zeroSolIndex)=0;
solIndexCount=length(solIndexNodes);
replacePattern(solIndexNodes)=1:solIndexCount;

for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        PHUTelem{indexPatch}(indexElem).solIndex=replacePattern(PHUTelem{indexPatch}(indexElem).nodesGlobal);
    end
end
end




