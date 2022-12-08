function [ m ] = ComputeMatrixM_vDirection_arbitraryDegree( m, solIndex,i,a,b,c,dRdu_patch1, dRdv_patch1,dRdu_patch2,  dRdv_patch2,p,q)
%Compute matrix m for v direction

solSet=length(solIndex)/3;
count=1;

basisIndex1=p;
basisIndex2=p+1;
basisIndex3=1;
basisIndex4=2;

for set=1:solSet
    
    index=solIndex(count:count+2);
    
    if index(1)~=0
        m(i,index(1))=b*dRdv_patch1(basisIndex1)+c*dRdu_patch1(basisIndex1);
    end
    
    if index(2)~=0
        m(i,index(2))=b*dRdv_patch1(basisIndex2)+c*dRdu_patch1(basisIndex2)+a*dRdu_patch2(basisIndex3);
    end
    
    if solIndex(3)~=0
        m(i,index(3))=a*dRdu_patch2(basisIndex4);
    end
  
    count=count+3;

    basisIndex1=basisIndex1+p+1;
    basisIndex2=basisIndex2+p+1;
    basisIndex3=basisIndex3+p+1;
    basisIndex4=basisIndex4+p+1;
 
end


end

