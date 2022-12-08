
close all
clear all
restoredefaultpath
addpaths

numPatches = 2;
L=1;
W=1;

numberElementsU=8;
numberElementsV=16;

GIFTmesh = init2DGeometryGIFTMP('rectangle',W,L,numPatches);

dimBasis = zeros(1, numPatches);
p=3;
q=3;

PHUTelem = cell(numPatches, 1);
for i=1:numPatches
    [PHUTelem{i}, dimBasis(i)] = initPHTmeshGenPeriodicBoundary( p,q,numberElementsU,numberElementsV );
    quadList{i} = 2:5;
end

for i=1:numPatches
    GIFTmesh{i}.c_net_rotated=GIFTmesh{i}.c_net;
end
patchBoundaries = {1, 2, 2, 4};

[PHUTelem, dimBasis, quadList,maxLevel ] = checkConforming_IGAPack( PHUTelem, dimBasis, patchBoundaries, p, q, quadList );
[PHUTelem, sizeBasis ] = zipConformingNew( PHUTelem, dimBasis, patchBoundaries, p, q);
[PHUTelem,sizeBasis ] = reassignNodesPeriodicBoundary( PHUTelem,numPatches,p,q,sizeBasis);


[ PHUTelem,solIndexCount ] = assignSolIndex_periodicBoundary( PHUTelem,patchBoundaries,p,q );

figure
plotPHTMeshMP(PHUTelem, GIFTmesh)
%pause
%
% figure
% plotPHTMesh_solIndex( PHUTelem,GIFTmesh,p)
% title('sol Index')
%
%
% figure
% plotPHTMesh_nodesGlobal( PHUTelem,GIFTmesh,p)
% title('nodes Global')
% %pause


%[PHUTelem,dimBasis,numType2Basis] =zipConformingC1PeriodicBoundary(PHUTelem,GIFTmesh,dimBasis, patchBoundaries, p, q,numConstraints);

disp('computing matrix m')
[m] = zipConformingC1_computeM(PHUTelem,GIFTmesh,patchBoundaries,solIndexCount,p,q);

disp('solving for coef sol')
[coefSol] = zipConformingC1_boundaryCondition_bilinearCahnHilliard(PHUTelem,m,p,q);

disp('modifying c')
[PHUTelem,sizeBasis,type2Basis] = zipConformingC1_modifyC(PHUTelem,GIFTmesh,patchBoundaries,sizeBasis,coefSol,p,q);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


% %%%%%%%%%%%%%%%%%%%%%%%
 disp('localizing type 2 basis function')
[PHUTelem] = localizeType2Basis2(PHUTelem,type2Basis);

% disp('finish localization')
% pause
% figure
% plotPHTMesh_nodesGlobal2( PHUTelem,GIFTmesh,type2basis,p)
% title('nodes Global')
%pause

% size(m)
% pause
%[PHUTelem,sizeBasis,numType2Basis] = zipConformingC1PeriodicBoundary(PHUTelem,GIFTmesh,sizeBasis, patchBoundaries, p, q,numConstraints );

%numType2Basis=length(type2Basis);
%testPlotBasisPhys_GIFTmesh( PHUTelem,GIFTmesh,p,q,sizeBasis,numType2Basis)
%[ PHUTelem] = localizeType2BasisPeriodicBoundary(PHUTelem);

%============================== CahnHilliard ==============================

n=sizeBasis;
maxR=0.05;
minR=-0.05;
initialSol=maxR + (minR-maxR).*rand(n,1);
initialSol=initialSol+0.63;
% %save('initialSolbillinearPatch.mat','initialSol')
% load initialSolbillinearPatch

%pause

plotSolPHUTCahnHilliard(PHUTelem, GIFTmesh,  p, q,initialSol)
az = 0;
el = 90;
view(az, el);
drawnow

saveas(gcf,['CahnHilliard2DC1_step',num2str(0),'.png'])

%================part on the cluster=====================%
C_n=initialSol;
Cdot_n=zeros(size(C_n));
timeStep=1e-11;
rhoInf=0.5; % rhoInfinity must be within [0,1];
alpha=3000;
theta=1.5;

tol=10e-3;
disp('Computing basis functions...')
tic

[ basisFun ] = saveBasisCahnHilliardBillinearPatch( PHUTelem,GIFTmesh, p, q );


toc

tolResidual=10e-5;
%start with 1st BE
disp(' Backward Eular')
BETimeStep=timeStep/2;
alphaM=1;
alphaF=1;
gamma=1;
C=C_n;
[C_n1,Cdot_n1] = InitializeGeneralizedAlphaMethod(C_n,Cdot_n,gamma);
[CalphaF,CalphaM] = GeneralizedAlphaMetohdRestartUpdate(C_n,Cdot_n,C_n1,Cdot_n1,alphaF,alphaM);

%[C_n1, Cdot_n1] = GeneralizedAlphaMethodSNESSolveBillinearPatch( PHUTelem,  GIFTmesh, basisFun,sizeBasis, p, q,theta,alpha,C_n,Cdot_n,C_n1,Cdot_n1,CalphaF,CalphaM,alphaM,alphaF,gamma,timeStep);
[C_n1, Cdot_n1] = GeneralizedAlphaMethodSNESSolveBillinearPatch_new( PHUTelem,  GIFTmesh, basisFun,sizeBasis, p, q,theta,alpha,C_n,Cdot_n,C_n1,Cdot_n1,CalphaF,CalphaM,alphaM,alphaF,gamma,timeStep);

C1=C_n1;
C_n=C_n1;
Cdot_n=Cdot_n1;
%2nd BE
[C_n1,Cdot_n1] = InitializeGeneralizedAlphaMethod(C_n,Cdot_n,gamma);
[CalphaF,CalphaM] = GeneralizedAlphaMetohdRestartUpdate(C_n,Cdot_n,C_n1,Cdot_n1,alphaF,alphaM);
%[C_n1, Cdot_n1] = GeneralizedAlphaMethodSNESSolveBillinearPatch( PHUTelem,  GIFTmesh, basisFun,sizeBasis, p, q,theta,alpha,C_n,Cdot_n,C_n1,Cdot_n1,CalphaF,CalphaM,alphaM,alphaF,gamma,timeStep);
[C_n1, Cdot_n1] = GeneralizedAlphaMethodSNESSolveBillinearPatch_new( PHUTelem,  GIFTmesh, basisFun,sizeBasis, p, q,theta,alpha,C_n,Cdot_n,C_n1,Cdot_n1,CalphaF,CalphaM,alphaM,alphaF,gamma,timeStep);


C2=C_n1;
%
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
%[C_n1, Cdot_n1] = GeneralizedAlphaMethodSNESSolveBillinearPatch( PHUTelem,  GIFTmesh, basisFun,sizeBasis, p, q,theta,alpha,C_n,Cdot_n,C_n1,Cdot_n1,CalphaF,CalphaM,alphaM,alphaF,gamma,timeStep);
[C_n1, Cdot_n1] = GeneralizedAlphaMethodSNESSolveBillinearPatch_new( PHUTelem,  GIFTmesh, basisFun,sizeBasis, p, q,theta,alpha,C_n,Cdot_n,C_n1,Cdot_n1,CalphaF,CalphaM,alphaM,alphaF,gamma,timeStep);

[error] = evaluateLTE( CBE,CLTE,C_n,C_n1,prevSol,step,timeStep,time,previousTime);
[timeStep,time,previousTime] = updateTimeStep(timeStep,time,error);
prevSol=C_n;
C_n=C_n1;
Cdot_n=Cdot_n1;

maxStep=10000;
for step=1:maxStep
    close all
    disp([num2str(step),' timeStep =',num2str(timeStep),', time =',num2str(time)])
    [C_n1,Cdot_n1] = InitializeGeneralizedAlphaMethod(C_n,Cdot_n,gamma);
    [CalphaF,CalphaM] = GeneralizedAlphaMetohdRestartUpdate(C_n,Cdot_n,C_n1,Cdot_n1,alphaF,alphaM);
   % [C_n1, Cdot_n1] = GeneralizedAlphaMethodSNESSolveBillinearPatch( PHUTelem,  GIFTmesh, basisFun,sizeBasis, p, q,theta,alpha,C_n,Cdot_n,C_n1,Cdot_n1,CalphaF,CalphaM,alphaM,alphaF,gamma,timeStep);
    [C_n1, Cdot_n1] = GeneralizedAlphaMethodSNESSolveBillinearPatch_new( PHUTelem,  GIFTmesh, basisFun,sizeBasis, p, q,theta,alpha,C_n,Cdot_n,C_n1,Cdot_n1,CalphaF,CalphaM,alphaM,alphaF,gamma,timeStep);

    
    [error] = evaluateLTE( CBE,CLTE,C_n,C_n1,prevSol,step,timeStep,time,previousTime);
    [timeStep,time,previousTime] = updateTimeStep(timeStep,time,error);
    prevSol=C_n;
    C_n=C_n1;
    Cdot_n=Cdot_n1;
    
    if round(step/10)==(step/10)
        plotSolPHUTCahnHilliard(PHUTelem, GIFTmesh,  p, q,C_n1)
        az = 0;
        el = 90;
        view(az, el);
        drawnow
        saveas(gcf,['CahnHilliard2DC1_step',num2str(step/10),'0.png'])
%        pause
    end
    
end
