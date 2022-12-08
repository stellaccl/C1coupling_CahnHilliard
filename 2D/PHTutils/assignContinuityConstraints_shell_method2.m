function [ m ,i] = assignContinuityConstraints_shell_method2(solIndex,basisIndex,m,i,dR_patchA,dR_patchB)
%Compute matrix m for v direction
%solSet=length(solIndex)/3;
%count=1;
%disp('in assign continuity constraints')
solSet=size(solIndex,2);

dRdu_patch1= dR_patchA(1,:);
dRdv_patch1= dR_patchA(2,:);
dRdw_patch1= dR_patchA(3,:);

dRdu_patch2= dR_patchB(1,:);
dRdv_patch2= dR_patchB(2,:);
dRdw_patch2= dR_patchB(3,:);


for set=1:solSet
    
    index=solIndex(:,set);
    
    basisIndex1=basisIndex(1,set);
    basisIndex2=basisIndex(2,set);
    basisIndex3=basisIndex(3,set);
    basisIndex4=basisIndex(4,set);
    
    if index(1)~=0
        m(i,index(1))=dRdu_patch1(basisIndex1);
        m(i+1,index(1))=dRdv_patch1(basisIndex1);
        m(i+2,index(1))=dRdw_patch1(basisIndex1);
    end
    
    if index(2)~=0
        
        m(i,index(2))=dRdu_patch1(basisIndex2)-dRdu_patch2(basisIndex3);
        m(i+1,index(2))=dRdv_patch1(basisIndex2)-dRdv_patch2(basisIndex3);
        m(i+2,index(2))=dRdw_patch1(basisIndex2)-dRdw_patch2(basisIndex3);
    end
    
    if index(3)~=0
        m(i,index(3))=-dRdu_patch2(basisIndex4);
        m(i+1,index(3))=-dRdv_patch2(basisIndex4);
        m(i+2,index(3))=-dRdw_patch2(basisIndex4);
    end
    %
end

i=i+2;
end

