restoredefaultpath
close all
clear

addpath('./utils/')
addpath('./nurbs_toolbox')



tic;

numPatches = 2;
W = 1;
L = 1;

GIFTmesh = init1DGeometryGIFTMP('straightLine');
%
p=2;
numElemU=3;

dimBasis = zeros(1, numPatches);
PHUTelem = cell(numPatches, 1);

for indexPatch=1:numPatches
    [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmesh_quadratic1D( p,numElemU);
end

figure
plotPHTMeshMP1D(PHUTelem, GIFTmesh)

%zipConformingC0
PHUTelem{1}(1).nodesGlobal=PHUTelem{1}(1).nodes;
PHUTelem{1}(2).nodesGlobal=PHUTelem{1}(2).nodes;
PHUTelem{1}(3).nodesGlobal=PHUTelem{1}(3).nodes;
PHUTelem{2}(1).nodesGlobal=[5,6,7];
PHUTelem{2}(2).nodesGlobal=[6,7,8];
PHUTelem{2}(3).nodesGlobal=[7,8,9];

%figure
for index=1:9
    plotBasis1D_GIFTmesh(PHUTelem,GIFTmesh,p,index)
end

solIndexCount=3;

PHUTelem{1}(1).solIndex=[0,0,0];
PHUTelem{1}(2).solIndex=[0,0,1];
PHUTelem{1}(3).solIndex=[0,1,2];

PHUTelem{2}(1).solIndex=[2,3,0];
PHUTelem{2}(2).solIndex=[3,0,0];
PHUTelem{2}(3).solIndex=[0,0,0];

%%% continuity condition
indexPatch=1;
indexElem=3;
uref=1;
[dRduPatch1] =computeBasisDerivatives1D( PHUTelem,GIFTmesh,indexPatch,indexElem,uref,p);

indexPatch=2;
indexElem=1;
uref=-1;
[dRduPatch2] =computeBasisDerivatives1D( PHUTelem,GIFTmesh,indexPatch,indexElem,uref,p);

m=[dRduPatch1(2) dRduPatch1(3)-dRduPatch2(1) -dRduPatch2(2)];
coefSol=null(m);

sizeBasis=9;
numType2Basis=size(coefSol,2);

%%%%%%%%%%% assignNewNodesGlobal %%%%%%%%%%%%%%
PHUTelem{1}(1).nodesGlobal=[1,2,3];
PHUTelem{1}(2).nodesGlobal=[2,3,7,8];
PHUTelem{1}(3).nodesGlobal=[3,7,8];
PHUTelem{2}(1).nodesGlobal=[7,8,4];
PHUTelem{2}(2).nodesGlobal=[7,4,5,8];
PHUTelem{2}(3).nodesGlobal=[4,5,6];
PHUTelem{1}(1).extraNodes=[];
PHUTelem{1}(2).extraNodes=[];
PHUTelem{1}(3).extraNodes=[];
PHUTelem{2}(1).extraNodes=[];
PHUTelem{2}(2).extraNodes=[];
PHUTelem{2}(3).extraNodes=[];

type2Basis=[7,8];

for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        PHUTelem{indexPatch}(indexElem).modifiedC=PHUTelem{indexPatch}(indexElem).C;
    end
end

%%%%%%%% modifyC  %%%%%%%%%
indexPatch=1;
indexElem=3;

solIndexInfo=zeros(solIndexCount,3);

solIndex=PHUTelem{indexPatch}(indexElem).solIndex;
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
        tempCoef=PHUTelem{patch}(elem).C(basis,:);
        mBezierCoef(:,indexColum)=tempCoef;
    end
end

mModifiedC=zeros((p+1),numType2Basis);
for indexBasis=1:numType2Basis
    sol=coefSol(:,indexBasis);
    mModifiedC(:,indexBasis)=mBezierCoef*sol;
end

PHUTelem{indexPatch}(indexElem).modifiedC(2,:)=mModifiedC(:,1)';
PHUTelem{indexPatch}(indexElem).modifiedC(3,:)=mModifiedC(:,2)';

%%%%%%%%%%%%%% modify c patch 1 elem 2 %%%%%%%%%%%%%%%%%%%%%
indexPatch=1;
indexElem=2;
%type2Loc=find(PHUTelem{indexPatch}(indexElem).solIndex~=0);

solIndexInfo=zeros(solIndexCount,3);

solIndex=PHUTelem{indexPatch}(indexElem).solIndex;
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
        tempCoef=PHUTelem{patch}(elem).C(basis,:);
        mBezierCoef(:,indexColum)=tempCoef;
    end
end

mModifiedC=zeros((p+1),numType2Basis);
for indexBasis=1:numType2Basis
    sol=coefSol(:,indexBasis);
    mModifiedC(:,indexBasis)=mBezierCoef*sol;
end

PHUTelem{indexPatch}(indexElem).modifiedC(3,:)=mModifiedC(:,1)';
PHUTelem{indexPatch}(indexElem).modifiedC(4,:)=mModifiedC(:,2)';

%%%%%%%%%%%%%%%%%%%%%

%%% modify c patch 2 elem 1
indexPatch=2;
indexElem=1;

solIndexInfo=zeros(solIndexCount,3);

solIndex=PHUTelem{indexPatch}(indexElem).solIndex;
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
        tempCoef=PHUTelem{patch}(elem).C(basis,:);
        mBezierCoef(:,indexColum)=tempCoef;
    end
end

mModifiedC=zeros((p+1),numType2Basis);
for indexBasis=1:numType2Basis
    sol=coefSol(:,indexBasis);
    mModifiedC(:,indexBasis)=mBezierCoef*sol;
end

PHUTelem{indexPatch}(indexElem).modifiedC(1,:)=mModifiedC(:,1)';
PHUTelem{indexPatch}(indexElem).modifiedC(2,:)=mModifiedC(:,2)';

%%%%%%%%%%%%%%%%%%%%%

%%% modify c patch 2 elem 2
indexPatch=2;
indexElem=2;
type2Loc=find(PHUTelem{indexPatch}(indexElem).solIndex~=0);

solIndexInfo=zeros(solIndexCount,3);

solIndex=PHUTelem{indexPatch}(indexElem).solIndex;
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
        tempCoef=PHUTelem{patch}(elem).C(basis,:);
        mBezierCoef(:,indexColum)=tempCoef;
    end
end

mModifiedC=zeros((p+1),numType2Basis);
for indexBasis=1:numType2Basis
    sol=coefSol(:,indexBasis);
    mModifiedC(:,indexBasis)=mBezierCoef*sol;
end

PHUTelem{indexPatch}(indexElem).modifiedC(1,:)=mModifiedC(:,1)';
PHUTelem{indexPatch}(indexElem).modifiedC(4,:)=mModifiedC(:,2)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %plot type2Basis
% %figure
figure
plotPHTMeshMP1D(PHUTelem, GIFTmesh)

[colors ] = colorGen([0, 0.5, 0],[1 1 0], 5);
color1=colors(1,:);
color2=colors(2,:);
color3=colors(3,:);
color4=colors(4,:);
color5=colors(5,:);
colorArray=[color1;color2;color3;color3;color2;color1;[1,0,0];[0,0,1]];

for index=1:8
    % plotBasis1D_modifiedC(PHUTelem,GIFTmesh,p,index)
    plotBasis1D_modifiedC_paper(PHUTelem,GIFTmesh,p,index,colorArray(index,:))
end
% 
% %type2basis function
hold on
line([0,1],[0 0],'LineWidth',2,'color','k');
hold on
line([0,0],[-0.1 0.1],'LineWidth',2,'color','k');
hold on
line([0+(1/6),0+(1/6)],[-0.05 0.05],'LineWidth',2,'color','k');
hold on
line([0+(2/6),0+(2/6)],[-0.05 0.05],'LineWidth',2,'color','k');
hold on
line([0+(3/6),0+(3/6)],[-0.1 0.1],'LineWidth',2,'color','k');
hold on
line([0+(4/6),0+(4/6)],[-0.05 0.05],'LineWidth',2,'color','k');
hold on
line([0+(5/6),0+(5/6)],[-0.05 0.05],'LineWidth',2,'color','k');
hold on
line([1,1],[-0.1 0.1],'LineWidth',2,'color','k');
%line([1,0],[-0.1 0.1],'LineWidth',2,'color','k');
Fsize=20;

adj=0;
text(0.02, 0.9+adj, '$${\omega_0}^{(1)}$$','interpreter','latex','Fontsize',Fsize)
hold on
text(0.1, 0.725+adj, '$${\omega_1}^{(1)}$$','interpreter','latex','Fontsize',Fsize)
hold on
text(0.25, 0.8+adj, '$${\omega_2}^{(1)}$$','interpreter','latex','Fontsize',Fsize)
hold on
text(0.43, 0.7+adj, '$${\psi_1}$$','interpreter','latex','Fontsize',Fsize)
hold on
text(0.6, 0.65+adj, '$${\psi_2}$$','interpreter','latex','Fontsize',Fsize)
hold on
text(0.75, 0.8+adj, '$${\omega_2}^{(2)}$$','interpreter','latex','Fontsize',Fsize)
hold on
text(0.85, 0.73+adj, '$${\omega_3}^{(2)}$$','interpreter','latex','Fontsize',Fsize)
hold on
text(0.895, 0.92+adj, '$${\omega_4}^{(2)}$$','interpreter','latex','Fontsize',Fsize)