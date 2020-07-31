function  [ newC,newCdot] = alphaStep(PHUTelem, GIFTmesh, sizeBasis, p, q, C,Cdot,timeStep,theta,alpha,alphaM,alphaF,gamma)
%referring to equations in Isogeometric Analysis of the CH ohase-field model
%use SNESSolve 

shift=alphaM/(alphaF*gamma*timeStep);
scale=1/alphaF;

C0=C;
Cdot0=Cdot;
Cnew=C0;
imax=100;
tol=1e-5; %for linear system

[newC,newCdot] = SNESSOlve( PHUTelem, GIFTmesh, sizeBasis, p, q,theta,alpha,C0,Cdot0,Cnew,alphaM,alphaF,gamma,timeStep,imax,shift,scale,tol);

end


