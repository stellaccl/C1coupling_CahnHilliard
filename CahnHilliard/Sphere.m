%Sphere.m
%sphere model

addpath('./exampleData')
addpath('./utils')
close all
clear

numPatches =6;
L=1;
W=1;

numberElementsU=12;
numberElementsV=12;
GIFTmesh = init2DGeometryGIFTMP('sphere',W,L,numPatches);

for indexPatch=1:numPatches
    GIFTmesh{indexPatch}.c_net_rotated=GIFTmesh{indexPatch}.c_net;
end

p=3;
q=3;

PHUTelem = cell(numPatches, 1);
dimBasis = zeros(1, numPatches);
for indexPatch=1:numPatches
    [PHUTelem{indexPatch}, dimBasis(indexPatch)] = initPHTmeshGenMF3( p,q,numberElementsU,numberElementsV );
    quadList{indexPatch} = 2:5;
end

patchBoundaries = {1,2,2,4;...
    2,3,2,2;...
    [1,3],4,[4,4],[4,2];...
    [1,2,3,4],5,[3,3,3,3],[1,2,3,4]
    [1,2,3,4],6,[1,1,1,1],[1,2,3,4]};
tic
disp('Zip conforming C0...')
[PHUTelem,sizeBasis] = zipConformingNew( PHUTelem, dimBasis, patchBoundaries, p, q);
toc
disp('Assign Solution Index...')
tic
[PHUTelem,solIndexCount,patchInfo] = assignSolIndex2D( PHUTelem,patchBoundaries,p,q);
toc

disp('Create boundary info...')
tic
[boundaryInfo,segInfo,subSegInfo,numSeg,patchNeighborInfo,patchEdgeInfo] = createBoundaryInfo3(PHUTelem, patchBoundaries );
toc
% figure
% plotPHTMesh_solIndex( PHUTelem,GIFTmesh)
disp('Zip conforming C1...')
tic
[PHUTelem,sizeBasis,numType2Basis,type2Basis] = zipConformingC1_3( PHUTelem,GIFTmesh,solIndexCount,sizeBasis, patchBoundaries,numPatches,numberElementsU,numberElementsV, p, q );
toc
disp('Localizing Type 2 Basis...')
tic
[PHUTelem] = localizeType2Basis2(PHUTelem,type2Basis);
toc

%============================== CahnHilliard ==============================
n=sizeBasis;
maxR=0.05;
minR=-0.05;
initialSol=maxR + (minR-maxR).*rand(n,1);
initialSol=initialSol+0.63;

plotSolPHUTCahnHilliardSphere(PHUTelem, GIFTmesh,  p, q,initialSol)
axis equal
saveas(gcf,['CahnHilliardSphere',num2str(numberElementsU),'_step',num2str(0),'.png'])

C_n=initialSol;
Cdot_n=zeros(size(C_n));
timeStep=1e-11;
rhoInf=0.5; % rhoInfinity must be within [0,1];
alpha=3000;
theta=1.5;

tol=10e-3;
disp('Computing basis functions...')
tic
[ basisFun ] = saveBasisCahnHilliard2( PHUTelem, p, q );
toc

tolResidual=10e-5;
%start with 1st BE
disp(' Backward Euler')
BETimeStep=timeStep/2;
alphaM=1;
alphaF=1;
gamma=1;
C=C_n;
[C_n1,Cdot_n1] = InitializeGeneralizedAlphaMethod(C_n,Cdot_n,gamma);
[CalphaF,CalphaM] = GeneralizedAlphaMetohdRestartUpdate(C_n,Cdot_n,C_n1,Cdot_n1,alphaF,alphaM);
[C_n1, Cdot_n1] = GeneralizedAlphaMethodSNESSolveSphere( PHUTelem, basisFun, sizeBasis, p, q,theta,alpha,C_n,Cdot_n,C_n1,Cdot_n1,CalphaF,CalphaM,alphaM,alphaF,gamma,timeStep);
C1=C_n1;
C_n=C_n1;
Cdot_n=Cdot_n1;

%2nd BE
[C_n1,Cdot_n1] = InitializeGeneralizedAlphaMethod(C_n,Cdot_n,gamma);
[CalphaF,CalphaM] = GeneralizedAlphaMetohdRestartUpdate(C_n,Cdot_n,C_n1,Cdot_n1,alphaF,alphaM);
[C_n2, Cdot_n2] = GeneralizedAlphaMethodSNESSolveSphere( PHUTelem, basisFun, sizeBasis, p, q,theta,alpha,C_n,Cdot_n,C_n1,Cdot_n1,CalphaF,CalphaM,alphaM,alphaF,gamma,timeStep);

C2=C_n1;
% Compute V0 ~ dX/dt at t0 with backward differences
Cdot_n=zeros(size(C));
[Cdot_n] = VecAXPY( Cdot_n,-3/timeStep,C);
[Cdot_n] = VecAXPY( Cdot_n,4/timeStep,C1);
[Cdot_n] = VecAXPY( Cdot_n,-1/timeStep,C2);

% Rough, lower-order estimate LTE of the initial step */
CLTE=zeros(size(C));
[CLTE] = VecAXPY( CLTE,2,C2);
[CLTE] = VecAXPY( CLTE,-4,C1);
[CLTE] = VecAXPY( CLTE,2,C);

CBE=C2;
% ================================== step 0 =================================
disp('Alpha step')
alphaM=0.5*((3-rhoInf)/(1+rhoInf));
alphaF=1/(1+rhoInf);
gamma=0.5+alphaM-alphaF;
time=0;
step=0;
previousTime=0;
prevSol=0;
C_n=C; %reassign tp initalSol
disp([num2str(step),' timeStep =',num2str(timeStep),', time =',num2str(time)])
[C_n1,Cdot_n1] = InitializeGeneralizedAlphaMethod(C_n,Cdot_n,gamma);
[CalphaF,CalphaM] = GeneralizedAlphaMetohdRestartUpdate(C_n,Cdot_n,C_n1,Cdot_n1,alphaF,alphaM);
[C_n1, Cdot_n1] = GeneralizedAlphaMethodSNESSolveSphere( PHUTelem, basisFun, sizeBasis, p, q,theta,alpha,C_n,Cdot_n,C_n1,Cdot_n1,CalphaF,CalphaM,alphaM,alphaF,gamma,timeStep);
[error] = evaluateLTE( CBE,CLTE,C_n,C_n1,prevSol,step,timeStep,time,previousTime);
[timeStep,time,previousTime] = updateTimeStep(timeStep,time,error);
prevSol=C_n;
C_n=C_n1;
Cdot_n=Cdot_n1;


maxStep=10000;
for step=1:maxStep
    disp([num2str(step),' timeStep =',num2str(timeStep),', time =',num2str(time)])
    [C_n1,Cdot_n1] = InitializeGeneralizedAlphaMethod(C_n,Cdot_n,gamma);
    [CalphaF,CalphaM] = GeneralizedAlphaMetohdRestartUpdate(C_n,Cdot_n,C_n1,Cdot_n1,alphaF,alphaM);
    [C_n1, Cdot_n1] = GeneralizedAlphaMethodSNESSolveSphere( PHUTelem, basisFun, sizeBasis, p, q,theta,alpha,C_n,Cdot_n,C_n1,Cdot_n1,CalphaF,CalphaM,alphaM,alphaF,gamma,timeStep);
    [error] = evaluateLTE( CBE,CLTE,C_n,C_n1,prevSol,step,timeStep,time,previousTime);
    [timeStep,time,previousTime] = updateTimeStep(timeStep,time,error);
    prevSol=C_n;
    C_n=C_n1;
    Cdot_n=Cdot_n1;
    
    if round(step/10)==(step/10)
        close all

        plotSolPHUTCahnHilliardSphere(PHUTelem, GIFTmesh,  p, q,C_n1)
        axis equal
        drawnow
        saveas(gcf,['CahnHilliardSphere',num2str(numberElementsU),'_step',num2str(step),'.png'])
      
    end
    
end
