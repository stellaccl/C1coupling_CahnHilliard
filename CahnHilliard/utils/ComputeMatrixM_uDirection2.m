function [ m ] = ComputeMatrixM_uDirection2( m, solIndex,basisIndex,i,a,b,c,dRdu_patch1, dRdv_patch1,dRdu_patch2,  dRdv_patch2,p,q)
%Compute matrix m for u direction
%use solIndex assigned at PHUTelem
%disp('in compute matrix m u direction2')
% basisIndex1=p^2;
% basisIndex2=(p^2)+p+1;
% basisIndex3=1;
% basisIndex4=p+2;

for set=1:p+1
    
    tempSolIndex= solIndex(:,set);
    tempBasisIndex=basisIndex(:,set);
    basisIndex1=tempBasisIndex(1);
    basisIndex2=tempBasisIndex(2);
    basisIndex3=tempBasisIndex(3);
    basisIndex4=tempBasisIndex(4);

    if tempSolIndex(1)~=0
        m(i,tempSolIndex(1))=b*dRdu_patch1(basisIndex1)+c*dRdv_patch1(basisIndex1);
    end
    
    if tempSolIndex(2)~=0
        m(i,tempSolIndex(2))=b*dRdu_patch1(basisIndex2)+c*dRdv_patch1(basisIndex2)+a*dRdv_patch2(basisIndex3);
    end
    
    if tempSolIndex(3)~=0
        m(i,tempSolIndex(3))=a*dRdv_patch2(basisIndex4);
    end
    
%     basisIndex1=basisIndex1+1;
%     basisIndex2=basisIndex2+1;
%     basisIndex3=basisIndex3+1;
%     basisIndex4=basisIndex4+1;
    
    
end


end

