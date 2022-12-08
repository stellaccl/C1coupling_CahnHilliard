function [ m ] = assignContinuityConstraints_u(solIndex,basisIndex,m,i,a,b,c,dRdu_patch1,dRdv_patch1,dRdu_patch2,dRdv_patch2,p,q)
%Compute matrix m for v direction
%solSet=length(solIndex)/3;
%count=1;

solSet=size(solIndex,2);

for set=1:solSet
    
    index=solIndex(:,set);
    
    basisIndex1=basisIndex(1,set);
    basisIndex2=basisIndex(2,set);
    basisIndex3=basisIndex(3,set);
    basisIndex4=basisIndex(4,set);
    
    if index(1)~=0
        m(i,index(1))=a*dRdu_patch1(basisIndex1)+b*dRdv_patch1(basisIndex1);
    end
    
    if index(2)~=0
        m(i,index(2))=a*dRdu_patch1(basisIndex2)+b*dRdv_patch1(basisIndex2)+c*dRdv_patch2(basisIndex3);
    end
    
    if index(3)~=0
        m(i,index(3))=c*dRdv_patch2(basisIndex4);
    end
    
    
end


end

