%OnePatch2D.m
%rectangle with periodic boundary conditions

addpath('./exampleData')
addpath('./utils')
close all
clear

numPatches = 1;
GIFTmesh = init2DGeometryGIFTMP('rectangle',1,1,numPatches);

tempDimBasis = zeros(1, numPatches);
dimBasis = zeros(1, numPatches);
p=3;
q=3;

numElemU=10;
numElemV=10;
PHUTelem = cell(numPatches, 1);
for i=1:numPatches
    [ PHUTelem{i}, dimBasis(i),tempDimBasis(i)] = initPHTmeshGenMF3( p,q, numElemU, numElemV );
    quadList{i} = 2:5;
end

%maxR=0.05;
%minR=-0.05;
%initialSol=maxR + (minR-maxR).*rand(n,1);
%initialSol=initialSol+0.63;
load('./exampleData/initialSol.mat')


[ PHUTelem{1}] =reassignNodesPR2( PHUTelem{1},numElemU,p);
dimBasis(1) = tempDimBasis(1);

plotSolPHUTCahnHilliard(PHUTelem, GIFTmesh,  p, q, initialSol )

C=initialSol;
Cdot=zeros(size(C));
timeStep=1e-11;
rhoInf=0.5; % rhoInfinity must be within [0,1];
alpha=3000;
theta=1.5;
%Backward Euler
alphaM=1;
alphaF=1;
gamma=1;
disp('alphaRestart ...')
[Cdot0,CLTE,CBE] = alphaRestart3(PHUTelem, GIFTmesh,dimBasis(1), p, q, C,Cdot,timeStep,theta,alpha,alphaM,alphaF,gamma);

alphaM=0.5*((3-rhoInf)/(1+rhoInf));
alphaF=1/(1+rhoInf);
gamma=0.5+alphaM-alphaF;
disp('alpha step ...')
time=0;
step=0;
previousTime=0;
prevSol=zeros(size(C));

disp([num2str(step),' timeStep =',num2str(timeStep),', time =',num2str(time)])
[newC,newCdot] = alphaStep(PHUTelem, GIFTmesh,dimBasis(1), p, q, C,Cdot0,timeStep,theta,alpha,alphaM,alphaF,gamma);
[error] = evaluateLTE(  CBE,CLTE,C,newC,prevSol,step,timeStep,time,previousTime);
[timeStep,time,previousTime] = updateTimeStep(timeStep,time,error);
prevSol=CBE;
maxStep=10000;
tol=10e-3;
for i=1:maxStep
    step=step+1
    C0=newC;
    Cdot0 = newCdot;
    disp([num2str(step),' timeStep =',num2str(timeStep),', time =',num2str(time)])
    [newC,newCdot] = alphaStep(PHUTelem, GIFTmesh,dimBasis(1), p, q,newC,newCdot,timeStep,theta,alpha,alphaM,alphaF,gamma);
    disp(['Norm of change in C: ', num2str(norm(newC-C0))])
    disp(['Norm of change in Cdot: ', num2str(norm(newCdot-Cdot0))])
    close all
    if round(i/10)==(i/10)
        plotSolPHUTCahnHilliard(PHUTelem, GIFTmesh,  p, q, newC)
        az = 0;
        el = 90;
        view(az, el);
        drawnow
        saveas(gcf,['CahnHilliard2D_step',num2str(i/10),'0.png'])
    end
    [error] = evaluateLTE( CBE,CLTE,C0,newC,prevSol,step,timeStep,time,previousTime );
    disp(['Error = ', num2str(error)])
    [timeStep,time,previousTime] = updateTimeStep(timeStep,time,error);
    prevSol=C0;
end







