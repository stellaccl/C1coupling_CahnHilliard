function [C_n1, Cdot_n1] = GeneralizedAlphaMethodSNESSolve2C1_2( PHUTelem, GIFTmesh, sizeBasis, p, q,theta,alpha,C_n,Cdot_n,C_n1,Cdot_n1,CalphaF,CalphaM,alphaM,alphaF,gamma,timeStep)
%

%disp('in SNESsolve ...')

repeat=1;
count=0;
tol=10e-5;
while repeat
    count=count+1
    shift=alphaM/(alphaF*gamma*timeStep);
    [residual] = assembleResidualCahnHilliardC1_2(PHUTelem, GIFTmesh, sizeBasis, p, q, CalphaF, CalphaM, theta, alpha);
    [tangent] = assembleTangentCahnHilliardC1_2(PHUTelem, GIFTmesh, sizeBasis, p, q, CalphaF, CalphaM, theta, alpha, shift);
    if norm(residual)<tol
        %update SNES result
        C_n1=C_n+(CalphaF-C_n)/alphaF;
        Cdot_n1=Cdot_n+(CalphaM-Cdot_n)/alphaM;
  
        break
    else

        %update CalphaF and CalphaM
        deltaCalphaF= tangent\(-residual);
%        cond(tangent)
%         tangent
%         norm(tangent)
        CalphaF=CalphaF+deltaCalphaF;
        CalphaM=(1-alphaM/gamma)*Cdot_n+(alphaM/(gamma*timeStep*alphaF))*(CalphaF-C_n);
        
    end
    
end


% disp('upadte after break')
% C_n1=C_n+(CalphaF-C_n)/alphaF;
% Cdot_n1=Cdot_n+(CalphaM-Cdot_n)/alphaM;

end

