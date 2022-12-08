function [C_n1, Cdot_n1] = GeneralizedAlphaMethodSNESSolveBillinearPatch_new( PHUTelem, GIFTmesh, basisFun, sizeBasis, p, q,theta,alpha,C_n,Cdot_n,C_n1,Cdot_n1,CalphaF,CalphaM,alphaM,alphaF,gamma,timeStep)
disp('in SNESsolve ...')
repeat=1;
count=0;
tol=10e-5;
oldNormResidual=1e3;

while repeat
    count=count+1
    shift=alphaM/(alphaF*gamma*timeStep);
    % [residual] = assembleResidualCahnHilliardC1(PHUTelem, GIFTmesh, sizeBasis, p, q, CalphaF, CalphaM, theta, alpha);
    % [residual] = assembleResidualCahnHilliardC1_check(PHUTelem, GIFTmesh, sizeBasis, p, q, CalphaF, CalphaM, theta, alpha);
    % [tangent] = assembleTangentCahnHilliardC1(PHUTelem, GIFTmesh, sizeBasis, p, q, CalphaF, CalphaM, theta, alpha, shift);
    [ tangent, residual ] = assembleTangentResidualCahnHilliard_BillinearPatch( PHUTelem, basisFun, sizeBasis, p, q, CalphaF, CalphaM, theta, alpha, shift);
    newNormResidual=norm(residual)
    
    %%%%%%%%%% new add 25012018%%%%%%%%%%%%%%%
    if newNormResidual>oldNormResidual || abs(newNormResidual-oldNormResidual)<10e-4
        disp('adjust time Step')
        timeStep=timeStep*2;
        if timeStep>10e-6
            disp('reset time step')
            timeStep=10e-8;
        end
        
        shift=alphaM/(alphaF*gamma*timeStep);
        [ tangent, residual ] = assembleTangentResidualCahnHilliard_BillinearPatch( PHUTelem, basisFun, sizeBasis, p, q, CalphaF, CalphaM, theta, alpha, shift);
        newNormResidual=norm(residual)
    end
    oldNormResidual= newNormResidual;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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

