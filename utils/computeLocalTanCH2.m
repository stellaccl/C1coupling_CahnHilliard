function localTangentTemp = computeLocalTanCH2(R,dR,laplacian,shift,dM,t1,t2,dLocalC,del2_c,M)
%compute the local Tangent matrix (vectorized)

vecA1 = dLocalC*dR;
vecB1 = t2*R' + dM*laplacian';
vecB2 = dM*del2_c*R'+M*laplacian';

localTangentTemp =  shift*(R*R') + t1*(dR'*dR) + vecA1'*vecB1 + laplacian*vecB2;

