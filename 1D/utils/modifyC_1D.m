function [PHTelem] = modifyC_1D(PHTelem,coefSol,solIndexCount,p )
% only for 2 patch domain

numType2Basis=size(coefSol,2);

for indexPatch=1:length(PHTelem)
    for indexElem=1:length(PHTelem{indexPatch})
        PHTelem{indexPatch}(indexElem).modifiedC=PHTelem{indexPatch}(indexElem).C;
    end
end


%%% patch 1 last element
indexPatch=1;
indexElem=length(PHTelem{indexPatch});

solIndexInfo=zeros(solIndexCount,3);

solIndex=PHTelem{indexPatch}(indexElem).solIndex;
for indexBasis=1:length(solIndex)
    if solIndex(indexBasis)~=0
        solIndexInfo(solIndex(indexBasis),:)=[indexPatch,indexElem,indexBasis];
    end
end

mBezierCoef=zeros((p+1),solIndexCount);
for indexColum=1:solIndexCount
    if all(solIndexInfo(indexColum,:))
        patch=solIndexInfo(indexColum,1);
        elem=solIndexInfo(indexColum,2);
        basis=solIndexInfo(indexColum,3);
        tempCoef=PHTelem{patch}(elem).C(basis,:);
        mBezierCoef(:,indexColum)=tempCoef;
    end
end

mModifiedC=zeros((p+1),numType2Basis);
for indexBasis=1:numType2Basis
    sol=coefSol(:,indexBasis);
    mModifiedC(:,indexBasis)=mBezierCoef*sol;
end

PHTelem{indexPatch}(indexElem).modifiedC(end-1,:)=mModifiedC(:,1)';
PHTelem{indexPatch}(indexElem).modifiedC(end,:)=mModifiedC(:,2)';

%%% patch 1 second last element
if p==2
    indexPatch=1;
    indexElem=length(PHTelem{indexPatch})-1;
    
    solIndexInfo=zeros(solIndexCount,3);
    
    solIndex=PHTelem{indexPatch}(indexElem).solIndex;
    for indexBasis=1:length(solIndex)
        if solIndex(indexBasis)~=0
            solIndexInfo(solIndex(indexBasis),:)=[indexPatch,indexElem,indexBasis];
        end
    end
    
    mBezierCoef=zeros((p+1),solIndexCount);
    for indexColum=1:solIndexCount
        if all(solIndexInfo(indexColum,:))
            patch=solIndexInfo(indexColum,1);
            elem=solIndexInfo(indexColum,2);
            basis=solIndexInfo(indexColum,3);
            tempCoef=PHTelem{patch}(elem).C(basis,:);
            mBezierCoef(:,indexColum)=tempCoef;
        end
    end
    
    mModifiedC=zeros((p+1),numType2Basis);
    for indexBasis=1:numType2Basis
        sol=coefSol(:,indexBasis);
        mModifiedC(:,indexBasis)=mBezierCoef*sol;
    end
    
    PHTelem{indexPatch}(indexElem).modifiedC(end,:)=mModifiedC(:,1)';
    PHTelem{indexPatch}(indexElem).modifiedC(end+1,:)=mModifiedC(:,2)';
end
%%%%%%%%%%%%%%%%%%%%% modify c patch 2 elem 1
indexPatch=2;
indexElem=1;

solIndexInfo=zeros(solIndexCount,3);

solIndex=PHTelem{indexPatch}(indexElem).solIndex;
for indexBasis=1:length(solIndex)
    if solIndex(indexBasis)~=0
        solIndexInfo(solIndex(indexBasis),:)=[indexPatch,indexElem,indexBasis];
    end
end

mBezierCoef=zeros((p+1),solIndexCount);
for indexColum=1:solIndexCount
    if all(solIndexInfo(indexColum,:))
        patch=solIndexInfo(indexColum,1);
        elem=solIndexInfo(indexColum,2);
        basis=solIndexInfo(indexColum,3);
        tempCoef=PHTelem{patch}(elem).C(basis,:);
        mBezierCoef(:,indexColum)=tempCoef;
    end
end

mModifiedC=zeros((p+1),numType2Basis);
for indexBasis=1:numType2Basis
    sol=coefSol(:,indexBasis);
    mModifiedC(:,indexBasis)=mBezierCoef*sol;
end

PHTelem{indexPatch}(indexElem).modifiedC(1,:)=mModifiedC(:,1)';
PHTelem{indexPatch}(indexElem).modifiedC(2,:)=mModifiedC(:,2)';

%%% modify c patch 2 elem 2
if p==2
    indexPatch=2;
    indexElem=2;
    %type2Loc=find(PHTelem{indexPatch}(indexElem).solIndex~=0);
    
    solIndexInfo=zeros(solIndexCount,3);
    
    solIndex=PHTelem{indexPatch}(indexElem).solIndex;
    for indexBasis=1:length(solIndex)
        if solIndex(indexBasis)~=0
            solIndexInfo(solIndex(indexBasis),:)=[indexPatch,indexElem,indexBasis];
        end
    end
    
    mBezierCoef=zeros((p+1),solIndexCount);
    for indexColum=1:solIndexCount
        if all(solIndexInfo(indexColum,:))
            patch=solIndexInfo(indexColum,1);
            elem=solIndexInfo(indexColum,2);
            basis=solIndexInfo(indexColum,3);
            tempCoef=PHTelem{patch}(elem).C(basis,:);
            mBezierCoef(:,indexColum)=tempCoef;
        end
    end
    
    mModifiedC=zeros((p+1),numType2Basis);
    for indexBasis=1:numType2Basis
        sol=coefSol(:,indexBasis);
        mModifiedC(:,indexBasis)=mBezierCoef*sol;
    end
    
    PHTelem{indexPatch}(indexElem).modifiedC(1,:)=mModifiedC(:,1)';
    PHTelem{indexPatch}(indexElem).modifiedC(end+1,:)=mModifiedC(:,2)';
end
end

