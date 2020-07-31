function [ PHUTelem ] = modifyC_arbitraryDegreePB3( PHUTelem,coefSol,numType2Basis,patchBoundaries,p,q)
%compute modifiedC
%use solIndex assogned at PHUTelem

numPatches = length(PHUTelem);

right_nodes1 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
right_nodes2 = (p+1):(p+1):(p+1)*(q+1);

left_nodes1 = 1:(p+1):(1+(p+1)*q);
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

%bezier coef index
coefIndex_patch1=sort([right_nodes1,right_nodes2]);
coefIndex_patch2=sort([left_nodes1,left_nodes2]);
numBoundaries = size(patchBoundaries,1);

% tempIndex1=[1,2];
% tempIndex2=[2,3];
%
% SolIndexPatch1=[1,2];
% SolIndexPatch2=[2,3];
%
% for i=1:p
%     tempIndex1=tempIndex1+3;
%     tempIndex2=tempIndex2+3;
%     SolIndexPatch1=[SolIndexPatch1,tempIndex1];
%     SolIndexPatch2=[SolIndexPatch2,tempIndex2];
% end


% SolIndexPatch1
% SolIndexPatch2

%assign index for basis functions which require bezier coef modification
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
    
    for indexElem = elementsA
        
        SolIndexPatch1=nonzeros(PHUTelem{1}(indexElem).solIndex);
        
        
        modifyBasisIndex=[right_nodes1,right_nodes2(1:end-1),right_nodes2(end):right_nodes2(end)+numType2Basis-(length(right_nodes2)*2)];
        modifyBasisIndex=sort(modifyBasisIndex);
        
        for i=1:numType2Basis
            
            localIndex =  modifyBasisIndex(i);
            PHUTelem{patchA}(indexElem).modifiedC(localIndex,:) = 0;
            PHUTelem{patchA}(indexElem).modifiedC(localIndex,coefIndex_patch1) = coefSol(SolIndexPatch1,i);
            
        end
        
        %         if indexElem == elementsA(end)
        %             for i=1:numType2Basis
        %                 localIndex =  modifyBasisIndex(i);
        %                 PHUTelem{patchA}(indexElem).modifiedC(localIndex,:) = 0;
        %                 SolIndexPatch1
        %                 tempCoefSol=coefSol(SolIndexPatch1,i)
        %                 tempCoefSol(end-p:end-p+1)=-tempCoefSol(end-p:end-p+1)
        %                 pause
        %                 PHUTelem{patchA}(indexElem).modifiedC(localIndex,coefIndex_patch1) =tempCoefSol;
        %             end
        %         else
        %             for i=1:numType2Basis
        %                 localIndex =  modifyBasisIndex(i);
        %                 PHUTelem{patchA}(indexElem).modifiedC(localIndex,:) = 0;
        %                 PHUTelem{patchA}(indexElem).modifiedC(localIndex,coefIndex_patch1) = coefSol(SolIndexPatch1,i);
        %             end
        %         end
        %
        
        
        
        
        %         SolIndexPatch1=SolIndexPatch1+(p*3)
        %         pause
        
    end
    
    for indexElem = elementsB
        SolIndexPatch2=nonzeros(PHUTelem{2}(indexElem).solIndex);
        modifyBasisIndex=[left_nodes1,left_nodes2,right_nodes2(end)+1:right_nodes2(end)+1+numType2Basis-(length(right_nodes2)*2+1)];
        modifyBasisIndex=sort(modifyBasisIndex);
        for i=1:numType2Basis
            localIndex =  modifyBasisIndex(i);
            PHUTelem{patchB}(indexElem).modifiedC(localIndex,:) = 0;
            PHUTelem{patchB}(indexElem).modifiedC(localIndex,coefIndex_patch2) = coefSol(SolIndexPatch2,i);
        end
        
        %         if indexElem == elementsB(end)
        %             for i=1:numType2Basis
        %                 localIndex =  modifyBasisIndex(i);
        %                 PHUTelem{patchB}(indexElem).modifiedC(localIndex,:) = 0;
        %                 SolIndexPatch2
        %                 tempCoefSol=coefSol(SolIndexPatch2,i);
        %                 tempCoefSol(end-p:end-p+1)=-tempCoefSol(end-p:end-p+1)
        %                 pause
        %                 PHUTelem{patchB}(indexElem).modifiedC(localIndex,coefIndex_patch2) =tempCoefSol;
        %             end
        %         else
        %             for i=1:numType2Basis
        %                 localIndex =  modifyBasisIndex(i);
        %                 PHUTelem{patchB}(indexElem).modifiedC(localIndex,:) = 0;
        %                 PHUTelem{patchB}(indexElem).modifiedC(localIndex,coefIndex_patch2) = coefSol(SolIndexPatch2,i);
        %             end
        %         end
        
        
        
        
        %         SolIndexPatch2=SolIndexPatch2+(p*3)
        %         pause
    end
    
end
end

