function localTangentTemp = computeLocalTanCH2(R,dR,laplacian,shift,dM,t1,t2,dLocalC,del2_c,M)
%compute the local Tangent matrix (vectorized)

vecA1 = dLocalC*dR;
vecB1 = t2*R' + dM*laplacian';
vecB2 = dM*del2_c*R'+M*laplacian';

localTangentTemp =  shift*(R*R') + t1*(dR'*dR) + vecA1'*vecB1 + laplacian*vecB2;

% 
% for iii=1:nument
%     Na = R(iii);
%     Na_x = dR(1,iii);
%     Na_y = dR(2,iii);
%     Na_z = dR(3,iii);
%     del2_Na=laplacian(iii);
%     for jjj=iii:nument
%         Nb = R(jjj);
%         Nb_x = dR(1,jjj);
%         Nb_y = dR(2,jjj);
%         Nb_z = dR(3,jjj);
%         del2_Nb=laplacian(jjj);
%         t3=t2*Nb+dM*del2_Nb;
%         
%         localTangentTemp(iii,jjj) =  shift*Na*Nb;
%         %localTangentTemp(iii,jjj) = Na*Nb;
%         localTangentTemp(iii,jjj) = localTangentTemp(iii,jjj) + (Na_x*Nb_x+Na_y*Nb_y+Na_z*Nb_z)*t1;
%         localTangentTemp(iii,jjj) = localTangentTemp(iii,jjj) + (Na_x*dLocalC(1)+Na_y*dLocalC(2)+Na_z*dLocalC(3))*t3;
%         localTangentTemp(iii,jjj) = localTangentTemp(iii,jjj) + del2_Na*(dM*del2_c*Nb+M*del2_Nb);
%         
%     end
% end