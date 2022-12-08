function  [localC, dLocalC, del2_c,localC_t]=computeConcentrationInfo2(tempC,tempCdot,R,dRdx,dRdy,d2Rdx,d2Rdy,d2Rdxdy,dxdxi,d2xdxi2)
%compute concentration C , its derivatives, and C_t at a given quadrature
%point (uref vref)



J=dxdxi';
G=J'*J;
dR=J*(inv(G))'*([dRdx'*tempC,dRdy'*tempC]');

dRdx_c=dRdx'*tempC;
dRdy_c=dRdy'*tempC;
d2Rdx_c=d2Rdx'*tempC;
d2Rdy_c=d2Rdy'*tempC;
d2Rdxdy_c=d2Rdxdy'*tempC;

[del2_c] = computeLaplaceBeltrami(G,dxdxi,d2xdxi2,dRdx_c,dRdy_c,d2Rdx_c,d2Rdy_c,d2Rdxdy_c);

localC=R'*tempC;
dcdu=dR(1);
dcdv=dR(2);
dcdw=dR(3);

localC_t=R'*tempCdot;
dLocalC=[dcdu,dcdv,dcdw];

end