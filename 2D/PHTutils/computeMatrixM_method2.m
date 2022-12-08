function [m,i] = computeMatrixM_method2(PHUTelem,GIFTmesh,patchInfo,m,i,p,q,numPts)

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
    
    %main boundary
    for ic=1:numPts
        
        [refPoint] = assignRefPoint(ic,uref,vref,edgeA,edgeB);
        
        urefA=refPoint.urefA;
        vrefA=refPoint.vrefA;
        urefB=refPoint.urefB;
        vrefB=refPoint.vrefB;
        
        i = i+1;
        
        [ dRdx_patchA, dRdy_patchA] = computeBasisDerivatives(PHUTelem,GIFTmesh,patchA,elemA(indexSeg),urefA,vrefA,p,q );
        [ dRdx_patchB, dRdy_patchB] = computeBasisDerivatives(PHUTelem,GIFTmesh,patchB,elemB(indexSeg),urefB,vrefB,p,q );

        testV1=ismember(edgeA,[2,4]) && ismember(edgeB,[2,4]);
        testV2=ismember(edgeA,[1,3]) && ismember(edgeB,[2,4]);
        
        testU1=ismember(edgeA,[1,3]) && ismember(edgeB,[1,3]);
        testU2=ismember(edgeA,[2,4]) && ismember(edgeB,[1,3]);
        
        if testV1 || testV2

            [a,b,c] = computeFunctionABC_GIFTmesh_vDirection(PHUTelem,GIFTmesh,patchA,elemA(indexSeg),patchB,elemB(indexSeg),refPoint);
            [m] = assignContinuityConstraints_v_method2(solIndex,basisIndex,m,i,a,b,c,dRdx_patchA,dRdy_patchA,dRdx_patchB,dRdy_patchB,p,q);
 
        elseif testU1 || testU2

            [a,b,c] = computeFunctionABC_GIFTmesh_uDirection(PHUTelem,GIFTmesh,patchA,elemA(indexSeg),patchB,elemB(indexSeg),refPoint);
            [m] = assignContinuityConstraints_u_method2(solIndex,basisIndex,m,i,a,b,c,dRdx_patchA,dRdy_patchA,dRdx_patchB,dRdy_patchB,p,q);

        else
            disp('error : undetermined case')
            pause
        end
        
    end
    
end




end

