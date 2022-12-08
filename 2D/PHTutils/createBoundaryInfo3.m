function [boundaryInfo,segInfo,subSegInfo,numSeg,patchNeighborInfo,patchEdgeInfo] = createBoundaryInfo3(PHUTelem, patchBoundaries )
%boundaryInfo{patchIndex} =[patchIndex,patchNeighborIndex,edgeOfPatchIndex,edgeOfPatchNeighbor,elemOfPatchIndex,elemOfPatchNeighbor]
%segInfo containts boundaryInfo all patches.
%with subBranch info
edgeU=[1,3];
edgeV=[2,4];

numPatches=length(PHUTelem);
segInfo=[];
tempBoundaryInfo=[];
numBoundaries = size(patchBoundaries,1);
for boundaryIndex = 1:numBoundaries
    
    patchAList = patchBoundaries{boundaryIndex,1};
    patchB = patchBoundaries{boundaryIndex,2};
    edgeAList = patchBoundaries{boundaryIndex,3};
    edgeBList = patchBoundaries{boundaryIndex,4};
    
    for indexPatch=1:length(patchAList)
        patchA = patchAList(indexPatch);
        edgeA = edgeAList(indexPatch);
        edgeB = edgeBList(indexPatch);
        
        [ elemA ] = sortEdgeElem( PHUTelem{patchA}, edgeA);
        [ elemB ] = sortEdgeElem( PHUTelem{patchB}, edgeB);
        
        elemA=cell2mat(elemA);
        elemB=cell2mat(elemB);
        
        for indexElem=1:length(elemA)
            segInfo=[ segInfo;patchA,patchB,edgeA,edgeB,elemA(indexElem),elemB(indexElem)];
            tempBoundaryInfo=[ tempBoundaryInfo;patchA,patchB,edgeA,edgeB,elemA(indexElem),elemB(indexElem)];
            tempBoundaryInfo=[ tempBoundaryInfo;patchB,patchA,edgeB,edgeA,elemB(indexElem),elemA(indexElem)];
        end
    end
end



%create patch info
patchNeighborInfo=cell(1,length(PHUTelem));
patchEdgeInfo=cell(1,length(PHUTelem));

%boundaryInfo=cell2mat(boundaryInfo);
for i=1:size(segInfo,1)
    patchNeighborInfo{segInfo(i,1)}=[patchNeighborInfo{segInfo(i,1)};segInfo(i,2)];
    patchNeighborInfo{segInfo(i,2)}=[patchNeighborInfo{segInfo(i,2)};segInfo(i,1)];
    patchEdgeInfo{segInfo(i,1)}=[patchEdgeInfo{segInfo(i,1)};segInfo(i,3)];
    patchEdgeInfo{segInfo(i,2)}=[patchEdgeInfo{segInfo(i,2)};segInfo(i,4)];
end

tempBoundaryInfo = sortrows(tempBoundaryInfo,1);
boundaryInfo=cell(1,numPatches);

for i=1:size(tempBoundaryInfo,1)
    boundaryInfo{tempBoundaryInfo(i,1)}=[boundaryInfo{tempBoundaryInfo(i,1)};tempBoundaryInfo(i,:)];
end

subBranchU=[];
subBranchV=[];

for indexPatch=1:length(PHUTelem)
    
    patchEdgeInfo{indexPatch}=unique(patchEdgeInfo{indexPatch});
    edge=patchEdgeInfo{indexPatch};
    for indexEdge=1:length(edge)
        [ elem ] = sortEdgeElem( PHUTelem{indexPatch}, edge(indexEdge));
        elem=cell2mat(elem);
        for indexElem=1:length(elem)-1
            
            if ismember(edge(indexEdge),edgeU)
                
                up=PHUTelem{indexPatch}(elem(indexElem)).neighbor_up;
                down=PHUTelem{indexPatch}(elem(indexElem)).neighbor_down;
                left=PHUTelem{indexPatch}(elem(indexElem)).neighbor_left;
                right=PHUTelem{indexPatch}(elem(indexElem)).neighbor_right;
                
                if ismember(elem(indexElem+1),up)
                    %disp('enter up')
                    edgeNeighbor=1;
                    edgeSelf=3;
                elseif ismember(elem(indexElem+1),down)
                    % disp('enter down')
                    edgeNeighbor=3;
                    edgeSelf=1;
                elseif ismember(elem(indexElem+1),right)
                    % disp('enter right')
                    edgeNeighbor=4;
                    edgeSelf=2;
                elseif ismember(elem(indexElem+1),left)
                    %disp('enter left')
                    edgeNeighbor=2;
                    edgeSelf=4;
                end
                %                 edgeNeighbor
                %                 edgeSelf
                %                 pause
                subBranchU=[subBranchU;indexPatch,indexPatch,edgeSelf,edgeNeighbor,elem(indexElem),elem(indexElem+1)];
                
            elseif ismember(edge(indexEdge),edgeV)
                
                up=PHUTelem{indexPatch}(elem(indexElem)).neighbor_up;
                down=PHUTelem{indexPatch}(elem(indexElem)).neighbor_down;
                left=PHUTelem{indexPatch}(elem(indexElem)).neighbor_left;
                right=PHUTelem{indexPatch}(elem(indexElem)).neighbor_right;
                
                if ismember(elem(indexElem+1),up)
                    %disp('enter up')
                    edgeNeighbor=1;
                    edgeSelf=3;
                elseif ismember(elem(indexElem+1),down)
                    % disp('enter down')
                    edgeNeighbor=3;
                    edgeSelf=1;
                elseif ismember(elem(indexElem+1),right)
                    % disp('enter right')
                    edgeNeighbor=4;
                    edgeSelf=2;
                elseif ismember(elem(indexElem+1),left)
                    %disp('enter left')
                    edgeNeighbor=2;
                    edgeSelf=4;
                end
                
                subBranchV=[subBranchV;indexPatch,indexPatch,edgeSelf,edgeNeighbor,elem(indexElem),elem(indexElem+1)];
            end
        end
    end
end

subSegInfo=struct;
subSegInfo.u=subBranchU;
subSegInfo.v=subBranchV;

numSeg=size(segInfo,1)+size(subBranchU,1)+size(subBranchV,1);


end

