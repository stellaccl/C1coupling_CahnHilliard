function [ m ,i] = assignConstraintsToM_shell( m,i,dR_patch1,dR_patch2,solIndex,basisIndex)
%Compute matrix m for v direction
%use information from paramMap
solSet=size(solIndex,1);

dRdu_patch1= dR_patch1(1,:);
dRdv_patch1= dR_patch1(2,:);
dRdw_patch1= dR_patch1(3,:);

dRdu_patch2= dR_patch2(1,:);
dRdv_patch2= dR_patch2(2,:);
dRdw_patch2= dR_patch2(3,:);

for set=1:solSet
    
    index=solIndex(set,:);
    basis=basisIndex(set,:);
    
    if index(1)~=0
      
        m(i,index(1))=dRdu_patch1(basis(1));
        m(i+1,index(1))=dRdv_patch1(basis(1));
        m(i+2,index(1))=dRdw_patch1(basis(1));
    end
    
    if index(2)~=0
        
        m(i,index(2))=dRdu_patch1(basis(2))-dRdu_patch2(basis(3));
        m(i+1,index(2))=dRdv_patch1(basis(2))-dRdv_patch2(basis(3));
        m(i+2,index(2))=dRdw_patch1(basis(2))-dRdw_patch2(basis(3));
    end
    
    if index(3)~=0
        
        m(i,index(3))=-dRdu_patch2(basis(4));
        m(i+1,index(3))=-dRdv_patch2(basis(4));
        m(i+2,index(3))=-dRdw_patch2(basis(4));
    end
    
end

i=i+2;
end

