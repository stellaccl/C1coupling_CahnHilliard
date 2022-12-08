function [C_alphaF,C_alphaM ] =GeneralizedAlphaMetohdRestartUpdate(C_n,Cdot_n,C_n1,Cdot_n1,alphaF,alphaM)

C_alphaF=C_n+alphaF*(C_n1-C_n);
C_alphaM=Cdot_n+alphaM*(Cdot_n1-Cdot_n);

end

