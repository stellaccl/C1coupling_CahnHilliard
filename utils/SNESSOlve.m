function [  newC,newCdot ] = SNESSOlve( PHUTelem, GIFTmesh, sizeBasis, p, q,theta,alpha,C0,Cdot0,newC,alphaM,alphaF,gamma,timeStep,imax,shift,scale,tol)
disp('in SNESsolve ...')
%newCdot = Cdot0;

for i=1:imax
    %i
    [C_alphaF,Cdot_alphaM,newCdot] = stageVector( C0,Cdot0,newC,alphaM,alphaF,gamma,timeStep);
    
    %assemble residual and tangent matrix
    [residual] = assembleResidualCahnHilliard( PHUTelem, GIFTmesh, sizeBasis, p, q,C_alphaF,Cdot_alphaM,theta,alpha);
    residual=residual*scale;
    
    if norm(residual)<tol
        %disp('norm(residual)<tol ...')
        break
    end

    [tangent] = assembleTangentCahnHilliard( PHUTelem, GIFTmesh, sizeBasis, p, q,C_alphaF,Cdot_alphaM,theta,alpha,shift);
    
    deltaCdot= tangent\(-residual);
    % residual
    % norm(residual)
    
    %update itarates using equation 27.1 and 27.2
    %newCdot=newCdot+deltaCdot;
    newC=newC+deltaCdot;
    
    
end

newC=C0+(newC-C0)/alphaF;
newCdot=Cdot0+(newCdot-Cdot0)/alphaM;

end

