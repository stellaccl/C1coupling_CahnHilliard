function  [localC, dLocalC, del2_c,localC_t]=computeConcentrationInfoC1_billinearPatch(tempC,tempCdot,R,dRdx,dRdy,d2Rdx,d2Rdy,d2Rdxdy,dxdxi,d2xdxi2, dxdxi2)
%compute concentration C , its derivatives, and C_t at a given quadrature
%point (uref vref)

localC=R'*tempC;
localC_t=R'*tempCdot;
dR = dxdxi\[dRdx'; dRdy'];
dcdu=dR(1,:)*tempC;
dcdv=dR(2,:)*tempC;
dLocalC=[dcdu,dcdv];
%ddR=[d2Rdx';d2Rdxdy';d2Rdy'];
ddR=[d2Rdx;d2Rdxdy;d2Rdy];

ddR = dxdxi2\(ddR - d2xdxi2*dR);


d2cdu=ddR(1,:)*tempC;
d2cdv=ddR(3,:)*tempC;
del2_c = d2cdu+d2cdv;
end