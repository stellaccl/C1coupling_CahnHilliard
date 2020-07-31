function [ m ] = ComputeMatrixM_vDirection2( m, solIndex,basisIndex,i,a,b,c,dRdu_patch1, dRdv_patch1,dRdu_patch2,  dRdv_patch2,p,q)
%Compute matrix m for v direction
%use sol index assignes to PHUTelem
%use pre-assigned basis index

% basisIndex1=p;
% basisIndex2=p+1;
% basisIndex3=1;
% basisIndex4=2;

% basisIndex1=2;
% basisIndex2=1;
% basisIndex3=p+1;
% basisIndex4=p;

%basisIndex
%pause

for set=1:p+1
    
    tempSolIndex= solIndex(:,set);
    tempBasisIndex=basisIndex(:,set);
    basisIndex1=tempBasisIndex(1);
    basisIndex2=tempBasisIndex(2);
    basisIndex3=tempBasisIndex(3);
    basisIndex4=tempBasisIndex(4);

    if tempSolIndex(1)~=0
        m(i,tempSolIndex(1))=b*dRdv_patch1(basisIndex1)+c*dRdu_patch1(basisIndex1);
    end
    
    if tempSolIndex(2)~=0
        m(i,tempSolIndex(2))=b*dRdv_patch1(basisIndex2)+c*dRdu_patch1(basisIndex2)+a*dRdu_patch2(basisIndex3);
    end
    
    if tempSolIndex(3)~=0
        m(i,tempSolIndex(3))=a*dRdu_patch2(basisIndex4);
    end
    
    %     basisIndex1=basisIndex1+p+1;
    %     basisIndex2=basisIndex2+p+1;
    %     basisIndex3=basisIndex3+p+1;
    %     basisIndex4=basisIndex4+p+1;
    
end


end

