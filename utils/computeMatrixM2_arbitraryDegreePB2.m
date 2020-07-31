function [m] = computeMatrixM2_arbitraryDegreePB2( PHUTelem,GIFTmesh,elementsA, elementsB, p,q ,numConstraints)
%Compute matrix m
%use solution index assigned in PHUTelem
numSeg=length(elementsA);
numColumns=numConstraints;

numPts=ceil(numColumns/numSeg);
fudge=0.1;
vref = linspace(-1+fudge,1-fudge,numPts);
uref = linspace(-1+fudge,1-fudge,numPts);
%number of colum = number of unknown coefs
%num of row in m -> num of equations needed to solve for the coefs
%m=zeros(numPts*(numSeg+(numSeg-1)*2),numColumns);
m=zeros(numPts*(numSeg+(numSeg-1)*2)+3,numColumns); %add 3 constraints
solIndex=1;
i=0;

right_nodes1 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
right_nodes2 = (p+1):(p+1):(p+1)*(q+1);

left_nodes1 = 1:(p+1):(1+(p+1)*q);
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;


%v direction
for indexSeg=1:numSeg
    
  %  solIndexList=solIndex:solIndex+((p+1)*3-1);
    solIndexList=[];
    for jj=1:p+1
       solIndexList=[solIndexList,PHUTelem{1}(elementsA(indexSeg)).solIndex([right_nodes1(jj),right_nodes2(jj)])];
       solIndexList=[solIndexList,PHUTelem{2}(elementsB(indexSeg)).solIndex([left_nodes2(jj)])];
    end
     
    for ic=1:numPts
        indexPatch1=1;
        indexElem1= elementsA(indexSeg);
        indexPatch2=2;
        indexElem2=elementsB(indexSeg);
        
        [a1,b1,c1] = computeFunctionABC_GIFTmesh_vDirection( PHUTelem,GIFTmesh,indexPatch1,indexElem1,indexPatch2,indexElem2,vref(ic));
        a1=-a1;   
        i = i+1;
        
        [ ~, dRdu_patch1,  dRdv_patch1] = tensorProductBernsteinPolynomial(1,vref(ic),p,q );
        [ ~, dRdu_patch2,  dRdv_patch2] = tensorProductBernsteinPolynomial(-1,vref(ic),p,q );
        
        [dRdu_patch1,dRdv_patch1,dRdu_patch2,dRdv_patch2] =scalingDerivatives2D(PHUTelem,indexPatch1,indexPatch2,indexElem1,indexElem2,dRdu_patch1,  dRdv_patch1, dRdu_patch2,  dRdv_patch2);

        [ m ] = ComputeMatrixM_vDirection_arbitraryDegree( m, solIndexList,i,a1,b1,c1,dRdu_patch1,  dRdv_patch1,dRdu_patch2,dRdv_patch2,p,q );
        
    end
    
    solIndex=solIndex+(p*3);
    
end

%u direction
%select bezier ordinate indices along the horizontal interface to ensure
%continuity
%assume numbering is as:
%...
%b7 b8 b9
%b4 b5 b6
%b1 b2 b3
%---------- --> down boundary
down_nodes1= 1:p+1;
down_nodes2=down_nodes1+(p+1);

up_nodes2=(p+1)*(q+1)-p:(p+1)*(q+1);
up_nodes1=up_nodes2-(p+1);

for indexSeg=1:numSeg-1
    
    solIndexList1=[];
    for jj=1:p+1
       solIndexList1=[solIndexList1,PHUTelem{1}(elementsA(indexSeg)).solIndex([up_nodes1(jj),up_nodes2(jj)])];
       solIndexList1=[solIndexList1,PHUTelem{1}(elementsA(indexSeg+1)).solIndex([down_nodes2(jj)])];
    end
  

    solIndexList2=[];
    for jj=1:p+1
       solIndexList2=[solIndexList2,PHUTelem{2}(elementsB(indexSeg)).solIndex([up_nodes1(jj),up_nodes2(jj)])];
       solIndexList2=[solIndexList2,PHUTelem{2}(elementsB(indexSeg+1)).solIndex([down_nodes2(jj)])];
    end

    
    for ic=1:numPts
        [ ~, dRdu_patch2,  dRdv_patch2] = tensorProductBernsteinPolynomial(uref(ic),-1,p,q );
        [ ~, dRdu_patch1,  dRdv_patch1] = tensorProductBernsteinPolynomial(uref(ic),1,p,q );
        
        indexPatch1=1;
        indexElem1=elementsA(indexSeg);
        indexPatch2=1;
        indexElem2=elementsA(indexSeg+1);
        
        [a1,b1,c1] = computeFunctionABC_GIFTmesh_uDirection( PHUTelem,GIFTmesh,indexPatch1,indexElem1,indexPatch2,indexElem2,uref(ic));
        
        a1=-a1;
        
        i = i+1;
        %  [ m ] = ComputeMatrixM_uDirection_quartic( m, solIndexList1,i,a1,b1,c1, dRdu_patch1, dRdv_patch1,dRdu_patch2, dRdv_patch2 );
        [ m ] = ComputeMatrixM_uDirection_arbitraryDegree( m, solIndexList1,i,a1,b1,c1, dRdu_patch1, dRdv_patch1,dRdu_patch2, dRdv_patch2,p,q );
        
        indexPatch1=2;
        indexElem1=elementsB(indexSeg);
        indexPatch2=2;
        indexElem2=elementsB(indexSeg+1);
        
        [a2,b2,c2] = computeFunctionABC_GIFTmesh_uDirection( PHUTelem,GIFTmesh,indexPatch1,indexElem1,indexPatch2,indexElem2,uref(ic));
        a2=-a2;
        i = i+1;
        % [ m ] = ComputeMatrixM_uDirection( m, solIndexList2,i,a2,b2,c2, dRdu_patch1, dRdv_patch1,dRdu_patch2, dRdv_patch2);
        [ m ] = ComputeMatrixM_uDirection_arbitraryDegree( m, solIndexList2,i,a2,b2,c2, dRdu_patch1, dRdv_patch1,dRdu_patch2, dRdv_patch2,p,q);
    end
    

  
end


end

