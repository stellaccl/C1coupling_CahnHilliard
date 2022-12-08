function [m,i] = computeMatrixM_shell_method2_multiElemGIFTmesh(PHUTelem,GIFTmesh,patchInfo,m,i,p,q,numPts)

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
        
        [coords,dxdxi_PatchA] = paramMapShell_multiElemGIFTmesh(elemA(indexSeg),GIFTmesh{patchA},urefA,vrefA, xminPatchA, yminPatchA, xmaxPatchA, ymaxPatchA);
        [coords,dxdxi_PatchB] = paramMapShell_multiElemGIFTmesh(elemB(indexSeg),GIFTmesh{patchB},urefB,vrefB, xminPatchB, yminPatchB, xmaxPatchB, ymaxPatchB);
        
        [ dRdx_patchA, dRdy_patchA] = computeBasisDerivatives(PHUTelem,GIFTmesh,patchA,elemA(indexSeg),urefA,vrefA,p,q );
        [ dRdx_patchB, dRdy_patchB] = computeBasisDerivatives(PHUTelem,GIFTmesh,patchB,elemB(indexSeg),urefB,vrefB,p,q );
        
        J=dxdxi_PatchA';
        G=J'*J;
        dR_patchA=J*(inv(G))'*([dRdx_patchA',dRdy_patchA']');
        
        J=dxdxi_PatchB';
        G=J'*J;
        dR_patchB=J*(inv(G))'*([dRdx_patchB',dRdy_patchB']');
        
        [m,i] = assignContinuityConstraints_shell_method2(solIndex,basisIndex,m,i,dR_patchA,dR_patchB);
        
    end
    
end




end

