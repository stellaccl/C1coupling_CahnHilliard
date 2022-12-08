function [m ] = computeMatrixM_shell(PHUTelem,GIFTmesh,solIndexCount, allSeg, numSeg,p, q)

right_nodes2 = (p+1)-1:(p+1):(p+1)*(q+1)-1;
right_nodes1 = (p+1):(p+1):(p+1)*(q+1);

left_nodes1 = 1:(p+1):(1+(p+1)*q);
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

down_nodes1=1:(p+1);
down_nodes2=down_nodes1+(p+1);

up_nodes1=(p+1)*(q+1)-p:(p+1)*(q+1);
up_nodes2=up_nodes1-(p+1);

numColumns=solIndexCount;
numPts=ceil(numColumns/numSeg);
if numPts<4
    numPts=4;
end
fudge=0;
vref = linspace(-1+fudge,1-fudge,numPts);
uref = linspace(-1+fudge,1-fudge,numPts);

m=zeros(numSeg*numPts*3,numColumns);

i=0;
for indexSeg=1:numSeg
    
    for ic=1:numPts
        
        indexPatch1=allSeg(indexSeg,1);
        indexPatch2=allSeg(indexSeg,2);
        indexEdge1=allSeg(indexSeg,3);
        indexEdge2=allSeg(indexSeg,4);
        indexElem1=allSeg(indexSeg,5);
        indexElem2=allSeg(indexSeg,6);
        
        xminPatch1 = PHUTelem{indexPatch1}(indexElem1).vertex(1);
        xmaxPatch1 = PHUTelem{indexPatch1}(indexElem1).vertex(3);
        yminPatch1 = PHUTelem{indexPatch1}(indexElem1).vertex(2);
        ymaxPatch1 = PHUTelem{indexPatch1}(indexElem1).vertex(4);
        
        switch  indexEdge1
            case 1
                uref_temp=uref(ic);
                vref_temp=-1;
                solIndexMidP1=PHUTelem{indexPatch1}(indexElem1).solIndex(down_nodes1);
                solIndexP1=PHUTelem{indexPatch1}(indexElem1).solIndex(down_nodes2);
                basisIndexMidP1=down_nodes1;
                basisIndexP1=down_nodes2;
            case 2
                uref_temp=1;
                vref_temp=vref(ic);
                solIndexMidP1=PHUTelem{indexPatch1}(indexElem1).solIndex(right_nodes1);
                solIndexP1=PHUTelem{indexPatch1}(indexElem1).solIndex(right_nodes2);
                basisIndexMidP1=right_nodes1;
                basisIndexP1=right_nodes2;
            case 3
                uref_temp=uref(ic);
                vref_temp=1;
                solIndexMidP1=PHUTelem{indexPatch1}(indexElem1).solIndex(up_nodes1);
                solIndexP1=PHUTelem{indexPatch1}(indexElem1).solIndex(up_nodes2);
                basisIndexMidP1=up_nodes1;
                basisIndexP1=up_nodes2;
            case 4
                uref_temp=-1;
                vref_temp=vref(ic);
                solIndexMidP1=PHUTelem{indexPatch1}(indexElem1).solIndex(left_nodes1);
                solIndexP1=PHUTelem{indexPatch1}(indexElem1).solIndex(left_nodes2);
                basisIndexMidP1=left_nodes1;
                basisIndexP1=left_nodes2;
        end
        
        % [~, dxdxi_Patch1]=paramMapSphere4(PHUTelem,indexPatch1,indexElem1,uref_temp,vref_temp);
        [coords,dxdxi_Patch1] = paramMapShell(GIFTmesh{indexPatch1}, uref_temp,vref_temp, xminPatch1, yminPatch1, xmaxPatch1, ymaxPatch1);
        [~,dRdu_patch1,dRdv_patch1] = tensorProductBernsteinPolynomial(uref_temp,vref_temp,p,q );
        
        dRdu_patch1 =squeeze(dRdu_patch1)*2/(xmaxPatch1-xminPatch1);
        dRdv_patch1 =squeeze(dRdv_patch1)*2/(ymaxPatch1-yminPatch1);
        
        J=dxdxi_Patch1';
        G=J'*J;
        dR_patch1=J*(inv(G))'*([dRdu_patch1,dRdv_patch1]');

%         quiver3(coords(1),coords(2),coords(3),norma(1,1),norma(1,2),norma(1,3))
%         disp('finish plot')
%        pause
        xminPatch2 = PHUTelem{indexPatch2}(indexElem2).vertex(1);
        xmaxPatch2 = PHUTelem{indexPatch2}(indexElem2).vertex(3);
        yminPatch2 = PHUTelem{indexPatch2}(indexElem2).vertex(2);
        ymaxPatch2 = PHUTelem{indexPatch2}(indexElem2).vertex(4);
        
        switch  indexEdge2
            case 1
                uref_temp=uref(ic);
                vref_temp=-1;
                solIndexMidP2=PHUTelem{indexPatch2}(indexElem2).solIndex(down_nodes1);
                solIndexP2=PHUTelem{indexPatch2}(indexElem2).solIndex(down_nodes2);
                basisIndexMidP2=down_nodes1;
                basisIndexP2=down_nodes2;
            case 2
                uref_temp=1;
                vref_temp=vref(ic);
                solIndexMidP2=PHUTelem{indexPatch2}(indexElem2).solIndex(right_nodes1);
                solIndexP2=PHUTelem{indexPatch2}(indexElem2).solIndex(right_nodes2);
                basisIndexMidP2=right_nodes1;
                basisIndexP2=right_nodes2;
            case 3
                uref_temp=uref(ic);
                vref_temp=1;
                solIndexMidP2=PHUTelem{indexPatch2}(indexElem2).solIndex(up_nodes1);
                solIndexP2=PHUTelem{indexPatch2}(indexElem2).solIndex(up_nodes2);
                basisIndexMidP2=up_nodes1;
                basisIndexP2=up_nodes2;
            case 4
                uref_temp=-1;
                vref_temp=vref(ic);
                solIndexMidP2=PHUTelem{indexPatch2}(indexElem2).solIndex(left_nodes1);
                solIndexP2=PHUTelem{indexPatch2}(indexElem2).solIndex(left_nodes2);
                basisIndexMidP2=left_nodes1;
                basisIndexP2=left_nodes2;
        end
        
        if solIndexMidP1 ~= solIndexMidP2
            disp('error diff sol index Mid')
            pause
        end
        [~, dxdxi_Patch2] = paramMapShell(GIFTmesh{indexPatch2}, uref_temp,vref_temp, xminPatch2, yminPatch2, xmaxPatch2, ymaxPatch2);
        
        %   [~, dxdxi_Patch2]=paramMapSphere4(PHUTelem,indexPatch2,indexElem2,uref_temp,vref_temp);
        [ ~,dRdu_patch2,  dRdv_patch2] = tensorProductBernsteinPolynomial(uref_temp,vref_temp,p,q );
        
        dRdu_patch2 =squeeze(dRdu_patch2)*2/(xmaxPatch2-xminPatch2);
        dRdv_patch2 =squeeze(dRdv_patch2)*2/(ymaxPatch2-yminPatch2);
        
        J=dxdxi_Patch2';
        G=J'*J;
        dR_patch2=J*(inv(G))'*([dRdu_patch2,dRdv_patch2]');
%         disp('compute matrix M sphere new 3')
%         size(dR_patch2)
%         pause
        solIndex=[solIndexP1',solIndexMidP1',solIndexP2'];
        
        basisIndex=[basisIndexP1',basisIndexMidP1',basisIndexMidP2',basisIndexP2'];
        
        i=i+1;
        [m,i] = assignConstraintsToM_shell(m,i,dR_patch1,dR_patch2,solIndex,basisIndex);
    end
    
end


end

