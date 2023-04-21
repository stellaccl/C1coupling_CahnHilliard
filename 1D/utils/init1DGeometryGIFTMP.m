function [  GIFTmesh] = init1DGeometryGIFTMP(object_type)

if strcmp(object_type, 'straightLine')
    
    GIFTmesh = cell(2,1);
    numberElementsU = 1;
    p = 1;
    
    patchIndex=1;
    crv = nrbline([0,0],[0.5,0]);
    
    knotU=crv.knots;
    coefs=crv.coefs;
    
    GIFTmesh{patchIndex} = genGIFTmesh1D( knotU, coefs, p, numberElementsU);
    
    
    patchIndex=2;
    
    crv = nrbline([0.5,0],[1,0]);
    
    knotU=crv.knots;
    coefs=crv.coefs;
    
    GIFTmesh{patchIndex} = genGIFTmesh1D( knotU, coefs, p, numberElementsU);
    
elseif  strcmp(object_type, 'straightLine1patch')
    
    GIFTmesh = cell(1,1);
    numberElementsU = 1;
    p = 1;
    
    patchIndex=1;
    crv = nrbline([0,0],[1,0]);
    
    knotU=crv.knots;
    coefs=crv.coefs;
    
    GIFTmesh{patchIndex} = genGIFTmesh1D( knotU, coefs, p, numberElementsU);
    
end
end

