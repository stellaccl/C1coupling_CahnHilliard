function [m,i] = computeMatrixM(PHUTelem,GIFTmesh,patchInfo,m,i,p,q,numPts)

fudge=0.1;
vref = linspace(-1+fudge,1-fudge,numPts);
uref = linspace(-1+fudge,1-fudge,numPts);

patchA=patchInfo.patchA;
patchB=patchInfo.patchB;
edgeA=patchInfo.edgeA;
edgeB=patchInfo.edgeB;
elemA=patchInfo.elemA;
elemB=patchInfo.elemB;

[nodesA,midNodesA,nodesB,midNodesB] = selectSolIndexNodesAndBasisNodes(edgeA,edgeB,p,q);

%dealing with segments in the main boundary
%disp('dealing with segments in main boundary')
for indexSeg=1:length(elemA)
    
    %check sol Index
    if PHUTelem{patchA}(elemA(indexSeg)).solIndex(midNodesA)~=PHUTelem{patchB}(elemB(indexSeg)).solIndex(midNodesB)
        disp('error: miss match sol index')
        pause
    end
    
    solIndex=[PHUTelem{patchA}(elemA(indexSeg)).solIndex(nodesA);...
        PHUTelem{patchA}(elemA(indexSeg)).solIndex(midNodesA);...
        PHUTelem{patchB}(elemB(indexSeg)).solIndex(nodesB)];
    
    basisIndex=[nodesA;midNodesA;midNodesB;nodesB];
    
    xminPatchA = PHUTelem{patchA}(elemA(indexSeg)).vertex(1);
    xmaxPatchA = PHUTelem{patchA}(elemA(indexSeg)).vertex(3);
    yminPatchA = PHUTelem{patchA}(elemA(indexSeg)).vertex(2);
    ymaxPatchA = PHUTelem{patchA}(elemA(indexSeg)).vertex(4);
    
    xminPatchB = PHUTelem{patchB}(elemB(indexSeg)).vertex(1);
    xmaxPatchB = PHUTelem{patchB}(elemB(indexSeg)).vertex(3);
    yminPatchB = PHUTelem{patchB}(elemB(indexSeg)).vertex(2);
    ymaxPatchB = PHUTelem{patchB}(elemB(indexSeg)).vertex(4);
    
    %main boundary
    for ic=1:numPts
        
        [refPoint] = assignRefPoint(ic,uref,vref,edgeA,edgeB);
        
        urefA=refPoint.urefA;
        vrefA=refPoint.vrefA;
        urefB=refPoint.urefB;
        vrefB=refPoint.vrefB;
        
        i = i+1;
        [ ~, dRdu_patchA, dRdv_patchA] = tensorProductBernsteinPolynomial(urefA,vrefA,p,q);
        [ ~, dRdu_patchB, dRdv_patchB] = tensorProductBernsteinPolynomial(urefB,vrefB,p,q);
        
        dRdu_patchA =squeeze(dRdu_patchA)*2/(xmaxPatchA-xminPatchA);
        dRdv_patchA =squeeze(dRdv_patchA)*2/(ymaxPatchA-yminPatchA);
        
        dRdu_patchB =squeeze(dRdu_patchB)*2/(xmaxPatchB-xminPatchB);
        dRdv_patchB =squeeze(dRdv_patchB)*2/(ymaxPatchB-yminPatchB);
        
        testV1=ismember(edgeA,[2,4]) && ismember(edgeB,[2,4]);
        testV2=ismember(edgeA,[1,3]) && ismember(edgeB,[2,4]);
        
        testU1=ismember(edgeA,[1,3]) && ismember(edgeB,[1,3]);
        testU2=ismember(edgeA,[2,4]) && ismember(edgeB,[1,3]);
        
        if testV1 || testV2
            [a,b,c] = computeFunctionABC_GIFTmesh_vDirection(PHUTelem,GIFTmesh,patchA,elemA(indexSeg),patchB,elemB(indexSeg),refPoint);
            [m] = assignContinuityConstraints_v(solIndex,basisIndex,m,i,a,b,c,dRdu_patchA,dRdv_patchA,dRdu_patchB,dRdv_patchB,p,q);
        elseif testU1 || testU2
            [a,b,c] = computeFunctionABC_GIFTmesh_uDirection(PHUTelem,GIFTmesh,patchA,elemA(indexSeg),patchB,elemB(indexSeg),refPoint);
            [m] = assignContinuityConstraints_u(solIndex,basisIndex,m,i,a,b,c,dRdu_patchA,dRdv_patchA,dRdu_patchB,dRdv_patchB,p,q);
        else
            disp('error : undetermine case')
            pause
        end
        
    end
    
end

%subBoundary
%use patch1 patch2 edge1 edge2 elem1 elem2
%subBoundary patchA
%disp('dealing with sub segments in patch A')
if length(elemA)>1
    %subBoundary for patch A
    
    %select sol index and basis index based on
    down=PHUTelem{patchA}(elemA(1)).neighbor_down;
    right=PHUTelem{patchA}(elemA(1)).neighbor_right;
    up=PHUTelem{patchA}(elemA(1)).neighbor_up;
    left=PHUTelem{patchA}(elemA(1)).neighbor_left;
    
    neighbor=elemA(2);
    
    if neighbor==down
        edge1=1;
        edge2=3;
    elseif neighbor==right
        edge1=2;
        edge2=4;
    elseif neighbor==up
        edge1=3;
        edge2=1;
    elseif neighbor==left
        edge1=4;
        edge2=2;
    else
        disp('error: undefine neghbor')
        pause
    end
    
    patch1=patchA;
    patch2=patchA;
    
    [nodes1,midNodes1,nodes2,midNodes2] = selectSolIndexNodesAndBasisNodes(edge1,edge2,p,q);
    
    for indexSubSeg=1:length(elemA)-1
        
        elem1=elemA(indexSubSeg);
        elem2=elemA(indexSubSeg+1);
        
        solIndex=[PHUTelem{patch1}(elem1).solIndex(nodes1);...
            PHUTelem{patch1}(elem1).solIndex(midNodes1);...
            PHUTelem{patch2}(elem2).solIndex(nodes2)];
        
        basisIndex=[nodes1;midNodes1;midNodes2;nodes2];
        
        xminPatch1 = PHUTelem{patch1}(elem1).vertex(1);
        xmaxPatch1 = PHUTelem{patch1}(elem1).vertex(3);
        yminPatch1 = PHUTelem{patch1}(elem1).vertex(2);
        ymaxPatch1 = PHUTelem{patch1}(elem1).vertex(4);
        
        xminPatch2 = PHUTelem{patch2}(elem2).vertex(1);
        xmaxPatch2 = PHUTelem{patch2}(elem2).vertex(3);
        yminPatch2 = PHUTelem{patch2}(elem2).vertex(2);
        ymaxPatch2 = PHUTelem{patch2}(elem2).vertex(4);
        
        for ic=1:numPts
            
            [refPoint] = assignRefPoint(ic,uref,vref,edge1,edge2);
            
            uref1=refPoint.urefA;
            vref1=refPoint.vrefA;
            uref2=refPoint.urefB;
            vref2=refPoint.vrefB;
            
            i = i+1;
            [ ~, dRdu_patch1, dRdv_patch1] = tensorProductBernsteinPolynomial(uref1,vref1,p,q);
            [ ~, dRdu_patch2, dRdv_patch2] = tensorProductBernsteinPolynomial(uref2,vref2,p,q);
            
            dRdu_patch1 =squeeze(dRdu_patch1)*2/(xmaxPatch1-xminPatch1);
            dRdv_patch1 =squeeze(dRdv_patch1)*2/(ymaxPatch1-yminPatch1);
            
            dRdu_patch2 =squeeze(dRdu_patch2)*2/(xmaxPatch2-xminPatch2);
            dRdv_patch2 =squeeze(dRdv_patch2)*2/(ymaxPatch2-yminPatch2);
            
            if ismember(edge1,[2,4]) && ismember(edge2,[2,4])
                
                [a,b,c] = computeFunctionABC_GIFTmesh_vDirection(PHUTelem,GIFTmesh,patch1,elem1,patch2,elem2,refPoint);
                [m] = assignContinuityConstraints_v(solIndex,basisIndex,m,i,a,b,c,dRdu_patch1,dRdv_patch1,dRdu_patch2,dRdv_patch2,p,q);
                
            elseif ismember(edge1,[1,3]) && ismember(edge2,[1,3])
                
                [a,b,c] = computeFunctionABC_GIFTmesh_uDirection(PHUTelem,GIFTmesh,patch1,elem1,patch2,elem2,refPoint);
                [m] = assignContinuityConstraints_u(solIndex,basisIndex,m,i,a,b,c,dRdu_patch1,dRdv_patch1,dRdu_patch2,dRdv_patch2,p,q);
                
            end
        end
    end
end

%subBoundary for patch B
%disp('dealing with sub segments in patchB')
if length(elemB)>1
    
    %select sol index and basis index based on
    down=PHUTelem{patchB}(elemB(1)).neighbor_down;
    right=PHUTelem{patchB}(elemB(1)).neighbor_right;
    up=PHUTelem{patchB}(elemB(1)).neighbor_up;
    left=PHUTelem{patchB}(elemB(1)).neighbor_left;
    
    neighbor=elemB(2);
    
    if neighbor==down
        edge1=1;
        edge2=3;
    elseif neighbor==right
        edge1=2;
        edge2=4;
    elseif neighbor==up
        edge1=3;
        edge2=1;
    elseif neighbor==left
        edge1=4;
        edge2=2;
    else
        disp('error: undefine neghbor')
        pause
    end
    
    patch1=patchB;
    patch2=patchB;
    
    [nodes1,midNodes1,nodes2,midNodes2] = selectSolIndexNodesAndBasisNodes(edge1,edge2,p,q);
    
    for indexSubSeg=1:length(elemB)-1
        
        elem1=elemB(indexSubSeg);
        elem2=elemB(indexSubSeg+1);
        
        solIndex=[PHUTelem{patch1}(elem1).solIndex(nodes1);...
            PHUTelem{patch1}(elem1).solIndex(midNodes1);...
            PHUTelem{patch2}(elem2).solIndex(nodes2)];
        
        basisIndex=[nodes1;midNodes1;midNodes2;nodes2];
        
        xminPatch1 = PHUTelem{patch1}(elem1).vertex(1);
        xmaxPatch1 = PHUTelem{patch1}(elem1).vertex(3);
        yminPatch1 = PHUTelem{patch1}(elem1).vertex(2);
        ymaxPatch1 = PHUTelem{patch1}(elem1).vertex(4);
        
        xminPatch2 = PHUTelem{patch2}(elem2).vertex(1);
        xmaxPatch2 = PHUTelem{patch2}(elem2).vertex(3);
        yminPatch2 = PHUTelem{patch2}(elem2).vertex(2);
        ymaxPatch2 = PHUTelem{patch2}(elem2).vertex(4);
        
        for ic=1:numPts
            
            [refPoint] = assignRefPoint(ic,uref,vref,edge1,edge2);
            
            uref1=refPoint.urefA;
            vref1=refPoint.vrefA;
            uref2=refPoint.urefB;
            vref2=refPoint.vrefB;
            
            i = i+1;
            [ ~, dRdu_patch1, dRdv_patch1] = tensorProductBernsteinPolynomial(uref1,vref1,p,q);
            [ ~, dRdu_patch2, dRdv_patch2] = tensorProductBernsteinPolynomial(uref2,vref2,p,q);
            
            dRdu_patch1 =squeeze(dRdu_patch1)*2/(xmaxPatch1-xminPatch1);
            dRdv_patch1 =squeeze(dRdv_patch1)*2/(ymaxPatch1-yminPatch1);
            
            dRdu_patch2 =squeeze(dRdu_patch2)*2/(xmaxPatch2-xminPatch2);
            dRdv_patch2 =squeeze(dRdv_patch2)*2/(ymaxPatch2-yminPatch2);
            
            if ismember(edge1,[2,4]) && ismember(edge2,[2,4])
                
                [a,b,c] = computeFunctionABC_GIFTmesh_vDirection(PHUTelem,GIFTmesh,patch1,elem1,patch2,elem2,refPoint);
                [m] = assignContinuityConstraints_v(solIndex,basisIndex,m,i,a,b,c,dRdu_patch1,dRdv_patch1,dRdu_patch2,dRdv_patch2,p,q);
                
            elseif ismember(edge1,[1,3]) && ismember(edge2,[1,3])
                
                [a,b,c] = computeFunctionABC_GIFTmesh_uDirection(PHUTelem,GIFTmesh,patch1,elem1,patch2,elem2,refPoint);
                [m] = assignContinuityConstraints_u(solIndex,basisIndex,m,i,a,b,c,dRdu_patch1,dRdv_patch1,dRdu_patch2,dRdv_patch2,p,q);
            end
        end
        
    end
    
    
end


end

