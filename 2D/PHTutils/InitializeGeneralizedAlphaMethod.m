function [ C_n1,Cdot_n1 ] =InitializeGeneralizedAlphaMethod(C_n,Cdot_n,gamma)
%initial C_n+1 amd Cdot_n+1
C_n1=C_n;
Cdot_n1=(gamma-1)/(gamma)*Cdot_n;

end

