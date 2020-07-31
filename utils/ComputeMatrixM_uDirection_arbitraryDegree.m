function [ m ] = ComputeMatrixM_uDirection_arbitraryDegree( m, solIndex,i,a,b,c,dRdu_patch1, dRdv_patch1,dRdu_patch2,  dRdv_patch2,p,q)
%Compute matrix m for u direction

numSolSet=length(solIndex)/3;

count=1;

basisIndex1=p^2;
basisIndex2=(p^2)+p+1;
basisIndex3=1;
basisIndex4=p+2;

for set=1:numSolSet
    
    index=solIndex(count:count+2);
    
    if index(1)~=0
        m(i,index(1))=b*dRdu_patch1(basisIndex1)+c*dRdv_patch1(basisIndex1);
    end
    
    if index(2)~=0
        m(i,index(2))=b*dRdu_patch1(basisIndex2)+c*dRdv_patch1(basisIndex2)+a*dRdv_patch2(basisIndex3);
    end
    
    if index(3)~=0
        m(i,index(3))=a*dRdv_patch2(basisIndex4);
    end
    
    count=count+3;
    
    basisIndex1=basisIndex1+1;
    basisIndex2=basisIndex2+1;
    basisIndex3=basisIndex3+1;
    basisIndex4=basisIndex4+1;
  
    
end


end

