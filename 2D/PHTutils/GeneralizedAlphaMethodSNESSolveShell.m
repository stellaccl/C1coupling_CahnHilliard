function [C_n1, Cdot_n1] = GeneralizedAlphaMethodSNESSolveShell( PHUTelem,basisFun, sizeBasis, p, q,theta,alpha,C_n,Cdot_n,C_n1,Cdot_n1,CalphaF,CalphaM,alphaM,alphaF,gamma,timeStep)
repeat=1;
count=0;
tol=10e-3;
while repeat
    count=count+1
    shift=alphaM/(alphaF*gamma*timeStep);
    % disp('Tangent and Residual..')
    %tic
    [tangent, residual] = assembleTangentResidualCahnHilliard( PHUTelem, basisFun, sizeBasis, p, q, CalphaF, CalphaM, theta, alpha, shift);
%     save('tangentFlat.mat','tangent')
%     save('residualFlat.mat','residual')
%     disp('finish save tangent residual')
%     pause
    
    disp(['norm_residual=',num2str(norm(residual))])
    %pause
    %   norm(residual)
    if norm(residual)<tol
        %update SNES result
        C_n1=C_n+(CalphaF-C_n)/alphaF;
        Cdot_n1=Cdot_n+(CalphaM-Cdot_n)/alphaM;
        
        break
    else
        
        %update CalphaF and CalphaM
        deltaCalphaF= tangent\(-residual);
        CalphaF=CalphaF+deltaCalphaF;
        CalphaM=(1-alphaM/gamma)*Cdot_n+(alphaM/(gamma*timeStep*alphaF))*(CalphaF-C_n);
        
    end
    
end

end

