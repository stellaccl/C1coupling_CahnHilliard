function [ GIFTmesh] = init2DGeometryGIFTMP(object_type,L,W,numPatches)
%creates a 2d GIFT mesh and associated control points, knotvectors, Bezier
%extraction operators, etc.

if strcmp(object_type, 'rectangle')
    
    GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %divide the patches along the x direction
    xVertices = linspace(0,L,numPatches+1);
    
    for patchIndex = 1:numPatches
        
        %set the dimensions of the patch
        patchMinX = xVertices(patchIndex);
        patchMaxX = xVertices(patchIndex+1);
        patchMinY = 0;
        patchMaxY = W;
        
        %initialize geometry on coarsest mesh
        coefs(1:3,1,1) = [patchMinX; patchMinY; 0];
        coefs(1:3,1,2) = [patchMinX; patchMaxY; 0];
        coefs(1:3,2,1) = [patchMaxX; patchMinY; 0];
        coefs(1:3,2,2) = [patchMaxX; patchMaxY; 0];
        coefs(4,1,1) = 1;
        coefs(4,1,2) = 1;
        coefs(4,2,1) = 1;
        coefs(4,2,2) = 1;
        
        knotU = [0 0 1 1];
        knotV = [0 0 1 1];
        
        GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
        
    end
    

    
elseif strcmp(object_type, 'sphere')
    
      GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    
    % ======================= patch 1 =========================
    patchIndex=1;
    
    coefs(1:3,1,1) = [0; 0; 0];
    coefs(1:3,1,2) = [1; 1; 0];
    coefs(1:3,2,1) = [4; 0; 0];
    coefs(1:3,2,2) = [3; 1; 0];
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    
    % ======================= patch 2 =========================
    patchIndex=2;
    
    coefs(1:3,1,1) = [4; 0; 0];
    coefs(1:3,1,2) = [3; 1; 0];
    coefs(1:3,2,1) = [4; 4; 0];
    coefs(1:3,2,2) = [3; 3; 0];
    
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    % ======================= patch 3 =========================
    patchIndex=3;
    
    coefs(1:3,1,1) = [0; 4; 0];
    coefs(1:3,1,2) = [1; 3; 0];
    coefs(1:3,2,1) = [4; 4; 0];
    coefs(1:3,2,2) = [3; 3; 0];
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    % ======================= patch 4 =========================
    patchIndex=4;
    
    coefs(1:3,1,1) = [0; 0; 0];
    coefs(1:3,1,2) = [1; 1; 0];
    coefs(1:3,2,1) = [0; 4; 0];
    coefs(1:3,2,2) = [1; 3; 0];
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    % ======================= patch 5 =========================
    
    patchIndex=5;
    
    coefs(1:3,1,1) = [1; 1; 0];
    coefs(1:3,1,2) = [1; 3; 0];
    coefs(1:3,2,1) = [3; 1; 0];
    coefs(1:3,2,2) = [3; 3; 0];
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    % ======================= patch 6 =========================
    patchIndex=6;
    
    coefs(1:3,1,1) = [8; 0; 0];
    coefs(1:3,1,2) = [8; 4; 0];
    coefs(1:3,2,1) = [4; 0; 0];
    coefs(1:3,2,2) = [4; 4; 0];
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %===================================================================================================
    
elseif strcmp(object_type, 'threePatches')
    
    vertices = [0,0; 2,0; 1,1; 2,1; 0,2; 1,2; 2,2];
    patches = [1,2,4,3; 3,4,7,6; 1,3,6,5];
    
    
    GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %divide the patches along the x direction
    %     xVertices = linspace(0,L,numPatches+1);
    
    for patchIndex = 1:numPatches
        
        SWcorner = vertices(patches(patchIndex,1),:);
        SEcorner = vertices(patches(patchIndex,2),:);
        NEcorner = vertices(patches(patchIndex,3),:);
        NWcorner = vertices(patches(patchIndex,4),:);
        
        coefs(1:3,1,1) = [SWcorner'; 0];
        coefs(1:3,1,2) = [NWcorner'; 0];
        coefs(1:3,2,1) = [SEcorner'; 0];
        coefs(1:3,2,2) = [NEcorner'; 0];
        
        coefs(4,1,1) = 1;
        coefs(4,1,2) = 1;
        coefs(4,2,1) = 1;
        coefs(4,2,2) = 1;
        
        knotU = [0 0 1 1];
        knotV = [0 0 1 1];
        
        GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
        
    end
    
end
