function  [ Cdot0,CLTE,C2] = alphaRestart3(PHUTelem, GIFTmesh, sizeBasis, p, q, C,Cdot,timeStep,theta,alpha,alphaM,alphaF,gamma)
%referring to equations in Isogeometric Analysis of the CH ohase-field model
%use SNESSolve


%first BE
timeStep=timeStep/2;
%stageTime=timeStep;
shift=alphaM/(alphaF*gamma*timeStep);
scale=1/alphaF;

C0=C;
Cdot0=Cdot;
newC=C0;
imax=10;
tol=1e-5;

%   th->stage_time = t + Alpha_f*dt;
%   th->shift_V = Alpha_m/(Alpha_f*Gamma*dt);
%   th->scale_F = 1/Alpha_f;
[C1,Cdot1] = SNESSOlve( PHUTelem, GIFTmesh, sizeBasis, p, q,theta,alpha,C0,Cdot0,newC,alphaM,alphaF,gamma,timeStep,imax,shift,scale,tol);

%2nd BE
%stageTime=stageTime+timeStep;
C0=C1; %ierr = VecCopy(X1,th->X0);CHKERRQ(ierr);
newC=C0; % ierr = VecCopy(th->X0,X2);CHKERRQ(ierr);
%Cdot0=zeros(size(C0));
Cdot0=((gamma-1)/gamma)*Cdot1;
[C2] = SNESSOlve( PHUTelem, GIFTmesh, sizeBasis, p, q,theta,alpha,C0,Cdot0,newC,alphaM,alphaF,gamma,timeStep,imax,shift,scale,tol);

% Compute V0 ~ dX/dt at t0 with backward differences
Cdot0=zeros(size(C0));
[Cdot0] = VecAXPY( Cdot0,-3/timeStep,C);
[Cdot0] = VecAXPY( Cdot0,4/timeStep,C1);
[Cdot0] = VecAXPY( Cdot0,-1/timeStep,C2);

% Rough, lower-order estimate LTE of the initial step */
CLTE=zeros(size(C0));
[CLTE] = VecAXPY( CLTE,2,C2);
[CLTE] = VecAXPY( CLTE,-4,C1);
[CLTE] = VecAXPY( CLTE,2,C);

end

