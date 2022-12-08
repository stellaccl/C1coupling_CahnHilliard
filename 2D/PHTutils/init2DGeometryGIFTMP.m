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
    
elseif strcmp(object_type, '2patch_curved')
    
    p = 2;
    q = 2;
    
    numberElementsU = 3;
    numberElementsV = 3;
    
    knotU = [0, 0, linspace(0,1,4), 1, 1];
    knotV = [0, 0, linspace(0,1,4), 1, 1];
    
    patchIndex=1;
    
    
    coefs(:,:,1) =[...
        
    -0.026261580765068  -0.192256303779966  -0.232323995542183  -0.157912567983780   0.071045670657458;...
    2.048077720502620   1.521473771627771   0.937630263092613   0.508333565640290   0.061865000289875;...
    0                   0                   0                   0                   0;...
    1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,2) =[...
    
0.393275652022517   0.417558639649381   0.296143701515062   0.304238030724017   0.271860713888199;...
2.278884246038668   1.550394617232756   1.008074560232799   0.692395721083571   0.368622552725388;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,3) =[...
    
1.015498405052568   0.987385957243808   1.000000000000000   0.943690038231428   0.911312721395610;...
2.432186010855915   1.586827874328400   1.000000000000000   0.684301391874616   0.303867919053751;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,4) =[...
    
1.587894001655664   1.500000000000000   1.500000000000000   1.477915766022430   1.461727107604521;...
2.219796399483548   1.500000000000000   1.000000000000000   0.595263770576117   0.174358651710478;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,5) =[...
    
2.039471396015333   1.879153614869196   1.849094030904296   1.899193337512463   1.979764176977613;...
1.937523020585575   1.506668983755331   1.025715640316920   0.534742435556875   0.036755055158250;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


surf1 = nrbmak(coefs, {knotU, knotV});
%nrbctrlplot(surf1)
surf1 = nrbpermute(surf1,[2,1]);
surf1.coefs=surf1.coefs(:,:,end:-1:1);

knotU = surf1.knots{1};
knotV = surf1.knots{2};
coefs = surf1.coefs;

GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

patchIndex=2;

coefs(:,:,1) =[...
    
2.039471396015333   1.879153614869196   1.849094030904296   1.899193337512463   1.979764176977613;...
1.937523020585575   1.506668983755331   1.025715640316920   0.534742435556875   0.036755055158250;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,2) =[...
    
2.450285710202310   2.460305571523943   2.470325432845577   2.460305571523943   2.470325432845577;...
1.957562743228842   1.466589538468797   1.045755362960187   0.745159523311180   0.374424654410738;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,3) =[...
    
2.971318498927255   3.000000000000000   2.991358221570522   2.971318498927255   2.901179469675820;...
1.967582604550476   1.500000000000000   1.055775224281821   0.654980771416478   0.334345209124204;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,4) =[...
    
3.501916651479254   3.500000000000000   3.500000000000000   3.412192397079132   3.322013645184430;...
1.960813774430296   1.500000000000000   1.000000000000000   0.540766161440768   0.240170321791761;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,5) =[...
    
4.088088606956585   4.108301324446237   4.051316920284611   3.923205324482445   3.659020704679810;...
1.901114532228672   1.446406840559980   0.978771283563889   0.503020656527832   0.104668881548526;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


surf2 = nrbmak(coefs, {knotU, knotV});
% hold on
% nrbctrlplot(srf2)

surf2 = nrbpermute(surf2,[2,1]);
surf2.coefs=surf2.coefs(:,:,end:-1:1);

knotU = surf2.knots{1};
knotV = surf2.knots{2};
coefs = surf2.coefs;

GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );

elseif strcmp(object_type, '2patch_curvedCommonBoundary')
    
    p = 2;
    q = 2;
    
    numberElementsU = 3;
    numberElementsV = 3;
    
    knotU = [0, 0, linspace(0,1,4), 1, 1];
    knotV = [0, 0, linspace(0,1,4), 1, 1];
    
    patchIndex=1;
    
    coefs(:,:,1) =[...
        
    0                   0                   0                   0                   0;...
    0   0.250000000000000   0.500000000000000   0.750000000000000   1.000000000000000;...
    0                   0                   0                   0                   0;...
    1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,2) =[...
    
0.125000000000000   0.125000000000000   0.125000000000000   0.125000000000000   0.125000000000000;...
0   0.250000000000000   0.500000000000000   0.750000000000000   1.000000000000000;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,3) =[...
    
0.250000000000000   0.250000000000000   0.250000000000000   0.250000000000000   0.250000000000000;...
0   0.250000000000000   0.500000000000000   0.750000000000000   1.000000000000000;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,4) =[...
    
0.375000000000000   0.357944898467390   0.302481153645892   0.312668372082494   0.375000000000000;...
0   0.250413148303262   0.501697869739434   0.752982591175605   1.000000000000000;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,5) =[...
    
0.500000000000000   0.476795780227741   0.371527856382858   0.386242727457950   0.500000000000000;...
0   0.252676974622507   0.503961696058678   0.755246417494850   1.000000000000000;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];

surf1 = nrbmak(coefs, {knotU, knotV});
%nrbctrlplot(surf1)
surf1 = nrbpermute(surf1,[2,1]);

knotU = surf1.knots{1};
knotV = surf1.knots{2};
coefs = surf1.coefs;

GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );

%%%%%%%%%%%%%%%%%%%%

patchIndex=2;


coefs(:,:,1) =[...
    
0.500000000000000   0.476795780227741   0.371527856382858   0.386242727457950   0.500000000000000;...
0   0.252676974622507   0.503961696058678   0.755246417494850   1.000000000000000;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,2) =[...
    
0.625000000000000   0.625000000000000   0.625000000000000   0.625000000000000   0.625000000000000;...
0   0.250000000000000   0.500000000000000   0.750000000000000   1.000000000000000;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,3) =[...
    
0.750000000000000   0.750000000000000   0.750000000000000   0.750000000000000   0.750000000000000;...
0   0.250000000000000   0.500000000000000   0.750000000000000   1.000000000000000;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,4) =[...
    
0.875000000000000   0.875000000000000   0.875000000000000   0.875000000000000   0.875000000000000;...
0   0.250000000000000   0.500000000000000   0.750000000000000   1.000000000000000;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,5) =[...
    
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000;...
0   0.250000000000000   0.500000000000000   0.750000000000000   1.000000000000000;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];

surf2 = nrbmak(coefs, {knotU, knotV});
% hold on
% nrbctrlplot(srf2)
surf2 = nrbpermute(surf2,[2,1]);

knotU = surf2.knots{1};
knotV = surf2.knots{2};
coefs = surf2.coefs;

GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );

elseif strcmp(object_type, '2patch_leftRight')
    
    GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %divide the patches along the x direction
    
    
    
    patchIndex=1;
    %set the dimensions of the patch
    patchMinX = 0;
    patchMaxX = 0.5;
    patchMinY = 0;
    patchMaxY = 1;
    
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    patchIndex=2;
    %set the dimensions of the patch
    patchMinX =0.5;
    patchMaxX =1;
    patchMinY =0;
    patchMaxY =1;
    tilt=0.2;
    %initialize geometry on coarsest mesh
    coefs(1:3,1,1) = [patchMinX; patchMinY; 0];
    coefs(1:3,1,2) = [patchMinX; patchMaxY; 0];
    coefs(1:3,2,1) = [patchMaxX; patchMinY; 0];
    coefs(1:3,2,2) = [patchMaxX; patchMaxY+tilt; 0];
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type, 'square_threePatchCurved')
    
    p = 2;
    q = 2;
    
    numberElementsU = 3;
    numberElementsV = 3;
    
    knotU = [0, 0, linspace(0,1,4), 1, 1];
    knotV = [0, 0, linspace(0,1,4), 1, 1];
    
    patchIndex=1;
    
    coefs(:,:,1) =[...
        
    0   0.113839285714286   0.219494047619048   0.343005952380952   0.500000000000000;...
    0   0.159970238095238   0.290922619047619   0.404017857142857   0.500000000000000;...
    0                   0                   0                   0                   0;...
    1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,2) =[...
    
0.250000000000000   0.375000000000000   0.500000000000000   0.618303571428571   0.622767857142857;...
0   0.125000000000000   0.250000000000000   0.350446428571429   0.465029761904762;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,3) =[...
    
0.500000000000000   0.625000000000000   0.749255952380953   0.747767857142857   0.752232142857143;...
0   0.125000000000000   0.225446428571429   0.328125000000000   0.418898809523810;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,4) =[...
    
0.750000000000000   0.875000000000000   0.871279761904762   0.874255952380953   0.874255952380953;...
0   0.125000000000000   0.220982142857143   0.343005952380953   0.426339285714286;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,5) =[...
    
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000;...
0   0.125000000000000   0.250000000000000   0.375000000000000   0.500000000000000;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];

surf1 = nrbmak(coefs, {knotU, knotV});
% hold on
% nrbctrlplot(srf2)
surf1 = nrbpermute(surf1,[2,1]);

knotU = surf1.knots{1};
knotV = surf1.knots{2};
coefs = surf1.coefs;

GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
%%%%%%%%%%%%%%%%%%%%%%%
patchIndex=2;

coefs(:,:,1) =[...
    
0.500000000000000   0.524553571428571   0.578125000000000   0.569196428571429   0.500000000000000;...
0.500000000000000   0.630208333333333   0.752232142857143   0.875744047619048   1.000000000000000;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,2) =[...
    
0.622767857142857   0.645089285714286   0.671875000000000   0.673363095238095   0.625000000000000;...
0.465029761904762   0.621279761904762   0.749255952380953   0.874255952380953   1.000000000000000;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,3) =[...
    
0.752232142857143   0.750000000000000   0.750000000000000   0.750000000000000   0.750000000000000;...
0.418898809523810   0.625000000000000   0.750000000000000   0.875000000000000   1.000000000000000;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,4) =[...
    
0.874255952380953   0.875000000000000   0.875000000000000   0.875000000000000   0.875000000000000;...
0.426339285714286   0.625000000000000   0.750000000000000   0.875000000000000   1.000000000000000;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,5) =[...
    
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000;...
0.500000000000000   0.625000000000000   0.750000000000000   0.875000000000000   1.000000000000000;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


surf2 = nrbmak(coefs, {knotU, knotV});
% hold on
% nrbctrlplot(srf2)
surf2 = nrbpermute(surf2,[2,1]);

knotU = surf2.knots{1};
knotV = surf2.knots{2};
coefs = surf2.coefs;

GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );

%%%%%%%%%%%%%%%%%%%%%%%
patchIndex=3;

coefs(:,:,1) =[...
    
0                   0                   0                   0                   0;...
0   0.250000000000000   0.500000000000000   0.750000000000000   1.000000000000000;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,2) =[...
    
0.113839285714286   0.125000000000000   0.125000000000000   0.125000000000000   0.125000000000000;...
0.159970238095238   0.375000000000000   0.625000000000000   0.875000000000000   1.000000000000000;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,3) =[...
    
0.219494047619048   0.250000000000000   0.250000000000000   0.250000000000000   0.250000000000000;...
0.290922619047619   0.500000000000000   0.750000000000000   0.875000000000000   1.000000000000000;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,4) =[...
    
0.343005952380952   0.375000000000000   0.375000000000000   0.375000000000000   0.375000000000000;...
0.404017857142857   0.625000000000000   0.750000000000000   0.875000000000000   1.000000000000000;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];


coefs(:,:,5) =[...
    
0.500000000000000   0.524553571428571   0.578125000000000   0.569196428571429   0.500000000000000;...
0.500000000000000   0.630208333333333   0.752232142857143   0.875744047619048   1.000000000000000;...
0                   0                   0                   0                   0;...
1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000];

surf3 = nrbmak(coefs, {knotU, knotV});
% hold on
% nrbctrlplot(srf2)
%surf2 = nrbpermute(surf2,[2,1]);

knotU = surf3.knots{1};
knotV = surf3.knots{2};
coefs = surf3.coefs;

GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );

elseif strcmp(object_type, 'squareWithHole')
    
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    w=sqrt(2)/2;
    
    % length of square
    L=8;
    
    %radius of circle (must be smaller than L/4)
    R=2;
    
    %r = 1;
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [R*cos(pi/4);R*sin(pi/4); 0];
    coefs(1:3,1,2) = [R*(cos(pi/4)+cos(3*pi/4))*w;R*(sin(pi/4)+sin(3*pi/4))*w; 0];
    coefs(1:3,1,3) = [R*cos(3*pi/4);R*sin(3*pi/4); 0];
    
    coefs(1:3,2,1) = [L/4; L/4; 0];
    coefs(1:3,2,2) = [0; L/2; 0];
    coefs(1:3,2,3) = [-L/4; L/4; 0];
    
    coefs(1:3,3,1) = [L/2; L/2; 0];
    coefs(1:3,3,2) = [0; L/2; 0];
    coefs(1:3,3,3) = [-L/2; L/2; 0];
    
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = w;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    srf1 = nrbmak(coefs, {knotU, knotV});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    indexPatch=1;
    % figure
    % nrbctrlplot(srf1)
    
    knotU = srf1.knots{1};
    knotV = srf1.knots{2};
    coefs = srf1.coefs;
    
    GIFTmesh{indexPatch} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    indexPatch=2;
    rotAngle = pi/2;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotz(rotAngle);
    T2 = vectrans(-pnt);
    
    srf2  = nrbtform(srf1 ,T*rotMatrix*T2);
    % hold on
    % nrbctrlplot(srf2)
    
    knotU = srf2.knots{1};
    knotV = srf2.knots{2};
    coefs = srf2.coefs;
    
    GIFTmesh{indexPatch} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    indexPatch=3;
    rotAngle = pi/2;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotz(rotAngle);
    T2 = vectrans(-pnt);
    
    srf3  = nrbtform(srf2 ,T*rotMatrix*T2);
    % hold on
    % nrbctrlplot(srf3)
    
    knotU = srf3.knots{1};
    knotV = srf3.knots{2};
    coefs = srf3.coefs;
    
    GIFTmesh{indexPatch} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    indexPatch=4;
    rotAngle = pi/2;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotz(rotAngle);
    T2 = vectrans(-pnt);
    
    srf4  = nrbtform(srf3 ,T*rotMatrix*T2);
    % hold on
    % nrbctrlplot(srf4)
    
    knotU = srf4.knots{1};
    knotV = srf4.knots{2};
    coefs = srf4.coefs;
    
    GIFTmesh{indexPatch} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type, 'squareWithHoleModified')
    
    %square with hole formed by 8 patches, four patch form the hold and four
    %patch at the corner
    
    
    w=sqrt(2)/2;
    
    % length of square
    L=8;
    
    %radius of circle (must be smaller than or equal L/4)
    R=sqrt(2)/2;
    %R=2;
    %r = 1;
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [R*cos(pi/4);R*sin(pi/4); 0];
    coefs(1:3,1,2) = [R*(cos(pi/4)+cos(3*pi/4))*w;R*(sin(pi/4)+sin(3*pi/4))*w; 0];
    coefs(1:3,1,3) = [R*cos(3*pi/4);R*sin(3*pi/4); 0];
    
    coefs(1:3,2,1) = [R*cos(pi/4); R*sin(pi/4)+(L/2-R*sin(pi/4))/2; 0];
    coefs(1:3,2,2) = [0; L/2; 0];
    coefs(1:3,2,3) = [R*cos(3*pi/4); R*sin(3*pi/4)+(L/2-R*sin(3*pi/4))/2; 0];
    
    
    coefs(1:3,3,1) = [R*cos(pi/4); L/2; 0];
    coefs(1:3,3,2) = [0; L/2; 0];
    coefs(1:3,3,3) = [R*cos(3*pi/4); L/2; 0];
    
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = w;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    srf1 = nrbmak(coefs, {knotU, knotV});
    %nrbctrlplot(srf1)
    %hold on
    
    indexPatch=1;
    
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    tempsrf1=srf1;
    
    
    tempsrf1 = nrbpermute(tempsrf1,[2,1]);
    
    knotU = tempsrf1.knots{1};
    knotV = tempsrf1.knots{2};
    coefs = tempsrf1.coefs;
    
    coefs=coefs(:,end:-1:1,:);
    
    GIFTmesh{indexPatch} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    patchIndex= 2;
    
    patchMinX = -L/2;
    patchMaxX = R*cos(3*pi/4);
    patchMinY = R*sin(3*pi/4);
    patchMaxY = L/2;
    
    %initialize geometry on coarsest mesh
    coefs=[];
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
    
    p=1;
    q=1;
    numberElementsU=1;
    numberElementsV=1;
    
    srf2 = nrbmak(coefs, {knotU, knotV});
    %     hold on
    %     nrbctrlplot(srf2)
    
    
    knotU = srf2.knots{1};
    knotV = srf2.knots{2};
    coefs = srf2.coefs;
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    patchIndex= 3;
    
    p=2;
    q=2;
    rotAngle = pi/2;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotz(rotAngle);
    T2 = vectrans(-pnt);
    
    srf3  = nrbtform(srf1 ,T*rotMatrix*T2);
    %     hold on
    %     nrbctrlplot(srf3)
    
    
    tempsrf3=srf3;
    
    
    %tempsrf3 = nrbpermute(tempsrf3,[2,1]);
    
    knotU = tempsrf3.knots{1};
    knotV = tempsrf3.knots{2};
    coefs = tempsrf3.coefs;
    
    coefs=coefs(:,end:-1:1,end:-1:1);
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    patchIndex= 4;
    
    tempCoord=srf3.coefs(:,1,3);
    
    patchMinX = -L/2;
    patchMaxX = tempCoord(1);
    patchMinY = -L/2;
    patchMaxY = tempCoord(2);
    
    %initialize geometry on coarsest mesh
    coefs=[];
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
    
    p=1;
    q=1;
    numberElementsU=1;
    numberElementsV=1;
    
    srf4 = nrbmak(coefs, {knotU, knotV});
    %     hold on
    %     nrbctrlplot(srf4)
    
    knotU = srf4.knots{1};
    knotV = srf4.knots{2};
    coefs = srf4.coefs;
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    patchIndex= 5;
    
    p=2;
    q=2;
    rotAngle = pi;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotz(rotAngle);
    T2 = vectrans(-pnt);
    
    srf5  = nrbtform(srf1 ,T*rotMatrix*T2);
    %     hold on
    %     nrbctrlplot(srf5)
    
    tempsrf5=srf5;
    
    tempsrf5 = nrbpermute(tempsrf5,[2,1]);
    
    knotU = tempsrf5.knots{1};
    knotV = tempsrf5.knots{2};
    coefs = tempsrf5.coefs;
    
    coefs=coefs(:,:,end:-1:1);
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    patchIndex= 6;
    
    tempCoord=srf5.coefs(:,1,3);
    
    patchMinX = tempCoord(1);
    patchMaxX = L/2;
    patchMinY = -L/2;
    patchMaxY = tempCoord(2);
    
    %initialize geometry on coarsest mesh
    coefs=[];
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
    
    p=1;
    q=1;
    numberElementsU=1;
    numberElementsV=1;
    
    srf6 = nrbmak(coefs, {knotU, knotV});
    %     hold on
    %     nrbctrlplot(srf6)
    %
    knotU = srf6.knots{1};
    knotV = srf6.knots{2};
    coefs = srf6.coefs;
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    patchIndex= 7;
    
    p=2;
    q=2;
    rotAngle = 3*pi/2;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotz(rotAngle);
    T2 = vectrans(-pnt);
    
    srf7 = nrbtform(srf1 ,T*rotMatrix*T2);
    %     hold on
    %     nrbctrlplot(srf7)
    
    knotU = srf7.knots{1};
    knotV = srf7.knots{2};
    coefs = srf7.coefs;
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    patchIndex= 8;
    
    tempCoord=srf7.coefs(:,1,3);
    
    patchMinX = tempCoord(1);
    patchMaxX = L/2;
    patchMinY = tempCoord(2);
    patchMaxY = L/2;
    
    %initialize geometry on coarsest mesh
    coefs=[];
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
    
    p=1;
    q=1;
    numberElementsU=1;
    numberElementsV=1;
    
    srf8 = nrbmak(coefs, {knotU, knotV});
    %     hold on
    %     nrbctrlplot(srf8)
    
    knotU = srf8.knots{1};
    knotV = srf8.knots{2};
    coefs = srf8.coefs;
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(object_type, 'squareWithTwoHoles')
    
    
    w=sqrt(2)/2;
    
    % length of square
    L=8;
    %length of retangle=L*2;
    
    %radius of circle (must be smaller than or equal L/4)
    R=sqrt(2)/2;
    
    %r = 1;
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [R*cos(pi/4);R*sin(pi/4); 0];
    coefs(1:3,1,2) = [R*(cos(pi/4)+cos(3*pi/4))*w;R*(sin(pi/4)+sin(3*pi/4))*w; 0];
    coefs(1:3,1,3) = [R*cos(3*pi/4);R*sin(3*pi/4); 0];
    
    coefs(1:3,2,1) = [R*cos(pi/4); R*sin(pi/4)+(L/2-R*sin(pi/4))/2; 0];
    coefs(1:3,2,2) = [0; L/2; 0];
    coefs(1:3,2,3) = [R*cos(3*pi/4); R*sin(3*pi/4)+(L/2-R*sin(3*pi/4))/2; 0];
    
    
    coefs(1:3,3,1) = [R*cos(pi/4); L/2; 0];
    coefs(1:3,3,2) = [0; L/2; 0];
    coefs(1:3,3,3) = [R*cos(3*pi/4); L/2; 0];
    
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = w;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    srf7 = nrbmak(coefs, {knotU, knotV});
    %nrbctrlplot(srf1)
    %hold on
    
    pnt = [-4,0,0];
    T  = vectrans(pnt);
    
    srf7  = nrbtform(srf7 ,T);
    % hold on
    nrbctrlplot(srf7)
    
    
    indexPatch=7;
    
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %     knotU = srf7.knots{1};
    %     knotV = srf7.knots{2};
    %     coefs = srf7.coefs;
    
    tempSrf7=srf7;
    tempSrf7 = nrbpermute(tempSrf7,[2,1]);
    
    tempSrf7.coefs=tempSrf7.coefs(:,end:-1:1,:);
    
    knotU = tempSrf7.knots{1};
    knotV = tempSrf7.knots{2};
    coefs = tempSrf7.coefs;
    
    GIFTmesh{indexPatch} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    patchIndex= 5;
    %
    p=2;
    q=2;
    rotAngle = pi/2;
    pnt = [-4,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotz(rotAngle);
    T2 = vectrans(-pnt);
    
    srf5  = nrbtform(srf7 ,T*rotMatrix*T2);
    hold on
    nrbctrlplot(srf5)
    
    tempSrf5=srf5;
    
    
    tempSrf5.coefs=tempSrf5.coefs(:,end:-1:1,end:-1:1);
    
    knotU = tempSrf5.knots{1};
    knotV = tempSrf5.knots{2};
    coefs = tempSrf5.coefs;
    
    %     knotU = srf5.knots{1};
    %     knotV = srf5.knots{2};
    %     coefs = srf5.coefs;
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    %
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    patchIndex= 3;
    %
    p=2;
    q=2;
    rotAngle = pi;
    pnt = [-4,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotz(rotAngle);
    T2 = vectrans(-pnt);
    
    srf3  = nrbtform(srf7 ,T*rotMatrix*T2);
    hold on
    nrbctrlplot(srf3)
    
    
    tempSrf3=srf3;
 
    tempSrf3.coefs=tempSrf3.coefs(:,end:-1:1,:);
    
    tempSrf3 = nrbpermute(tempSrf3,[2,1]);
     
    
    knotU = tempSrf3.knots{1};
    knotV = tempSrf3.knots{2};
    coefs = tempSrf3.coefs;
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    %
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    patchIndex= 6;
    
    tempCoord=srf7.coefs(:,1,3);
    
    patchMinX = -L;
    patchMaxX = tempCoord(1);
    patchMinY = tempCoord(2);
    patchMaxY = L/2;
    
    %initialize geometry on coarsest mesh
    coefs=[];
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
    
    p=1;
    q=1;
    numberElementsU=1;
    numberElementsV=1;
    
    srf6 = nrbmak(coefs, {knotU, knotV});
    hold on
    nrbctrlplot(srf6)
    
    
    knotU = srf6.knots{1};
    knotV = srf6.knots{2};
    coefs = srf6.coefs;
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    patchIndex= 4;
    
    tempCoord=srf5.coefs(:,1,3);
    
    patchMinX = -L;
    patchMaxX = tempCoord(1);
    patchMinY = tempCoord(2);
    patchMaxY = -L/2;
    
    %initialize geometry on coarsest mesh
    coefs=[];
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
    
    p=1;
    q=1;
    numberElementsU=1;
    numberElementsV=1;
    
    srf4 = nrbmak(coefs, {knotU, knotV});
    hold on
    nrbctrlplot(srf4)
    
    tempSrf4=srf4;
 
    tempSrf4.coefs=tempSrf4.coefs(:,:,end:-1:1);

    knotU = tempSrf4.knots{1};
    knotV = tempSrf4.knots{2};
    coefs = tempSrf4.coefs;

    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    indexPatch=9;
    
    L=8;
    %length of retangle=L*2;
    
    %radius of circle (must be smaller than or equal L/4)
    %R=sqrt(2)/2;
    R=1.2;
    %r = 1;
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [R*cos(pi/4);R*sin(pi/4); 0];
    coefs(1:3,1,2) = [R*(cos(pi/4)+cos(3*pi/4))*w;R*(sin(pi/4)+sin(3*pi/4))*w; 0];
    coefs(1:3,1,3) = [R*cos(3*pi/4);R*sin(3*pi/4); 0];
    
    coefs(1:3,2,1) = [R*cos(pi/4); R*sin(pi/4)+(L/2-R*sin(pi/4))/2; 0];
    coefs(1:3,2,2) = [0; L/2; 0];
    coefs(1:3,2,3) = [R*cos(3*pi/4); R*sin(3*pi/4)+(L/2-R*sin(3*pi/4))/2; 0];
    
    
    coefs(1:3,3,1) = [R*cos(pi/4); L/2; 0];
    coefs(1:3,3,2) = [0; L/2; 0];
    coefs(1:3,3,3) = [R*cos(3*pi/4); L/2; 0];
    
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = w;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    srf9 = nrbmak(coefs, {knotU, knotV});
    
    nrbctrlplot(srf9)
    hold on
    
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
%     knotU = srf9.knots{1};
%     knotV = srf9.knots{2};
%     coefs = srf9.coefs;
    

    tempSrf9=srf9;
 
    tempSrf9 = nrbpermute(tempSrf9,[2,1]);
    tempSrf9.coefs=tempSrf9.coefs(:,end:-1:1,:);

    knotU = tempSrf9.knots{1};
    knotV = tempSrf9.knots{2};
    coefs = tempSrf9.coefs;


    GIFTmesh{indexPatch} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    patchIndex= 1;
    
    p=2;
    q=2;
    rotAngle = 3*pi/2;
    pnt = [-4,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotz(rotAngle);
    T2 = vectrans(-pnt);
    
    tempSrf1  = nrbtform(srf7 ,T*rotMatrix*T2);
    % hold on
    % nrbctrlplot(tempSrf1)
    
    rotAngle = pi/2;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotz(rotAngle);
    T2 = vectrans(-pnt);
    
    tempSrf2  = nrbtform(srf9 ,T*rotMatrix*T2);
    % hold on
    % nrbctrlplot(tempSrf2)
    
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,:) =tempSrf1.coefs(1:3,1,:);
    coefs(1:3,3,:) =tempSrf2.coefs(1:3,1,end:-1:1);
    
    P1=(coefs(1:3,1,1)+coefs(1:3,3,1))/2;
    P2=(coefs(1:3,1,2)/w+coefs(1:3,3,2)/w)/2;
    P3=(coefs(1:3,1,3)+coefs(1:3,3,3))/2;
    
    %check=coefs(1:3,1,2)*w;
    
    coefs(1:3,2,1) = P1;
    coefs(1:3,2,2) = P2;
    coefs(1:3,2,3) = P3;
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = w;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = w;
    coefs(4,3,3) = 1;
    
    srf1 = nrbmak(coefs, {knotU, knotV});
    nrbctrlplot(srf1)
    
    knotU = srf1.knots{1};
    knotV = srf1.knots{2};
    coefs = srf1.coefs;
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    patchIndex= 13;
    %
    p=2;
    q=2;
    rotAngle = pi;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotz(rotAngle);
    T2 = vectrans(-pnt);
    
    srf13  = nrbtform(srf9 ,T*rotMatrix*T2);
    hold on
    nrbctrlplot(srf13)
    
%     knotU = srf13.knots{1};
%     knotV = srf13.knots{2};
%     coefs = srf13.coefs;
    

     tempSrf13=srf13;
 
    tempSrf13 = nrbpermute(tempSrf13,[2,1]);
    tempSrf13.coefs=tempSrf13.coefs(:,:,end:-1:1);

    knotU = tempSrf13.knots{1};
    knotV = tempSrf13.knots{2};
    coefs = tempSrf13.coefs;

    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    patchIndex= 11;
    %
    p=2;
    q=2;
    rotAngle = 3*pi/2;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotz(rotAngle);
    T2 = vectrans(-pnt);
    
    srf11  = nrbtform(srf9 ,T*rotMatrix*T2);
    hold on
    nrbctrlplot(srf11)
    
    knotU = srf11.knots{1};
    knotV = srf11.knots{2};
    coefs = srf11.coefs;
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    patchIndex= 10;
    
    tempCoord=srf9.coefs(:,1,1);
    
    patchMinX = tempCoord(1);
    patchMaxX = L/2;
    patchMinY = tempCoord(1);
    patchMaxY = L/2;
    
    %initialize geometry on coarsest mesh
    coefs=[];
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
    
    p=1;
    q=1;
    numberElementsU=1;
    numberElementsV=1;
    
    srf10 = nrbmak(coefs, {knotU, knotV});
    hold on
    nrbctrlplot(srf10)
    
    
    knotU = srf10.knots{1};
    knotV = srf10.knots{2};
    coefs = srf10.coefs;
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    patchIndex= 12;
    
    tempCoord=srf13.coefs(:,1,3);
    
    patchMinX = tempCoord(1);
    patchMaxX = L/2;
    patchMinY =-L/2;
    patchMaxY =tempCoord(2);
    
    %initialize geometry on coarsest mesh
    coefs=[];
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
    
    p=1;
    q=1;
    numberElementsU=1;
    numberElementsV=1;
    
    srf12 = nrbmak(coefs, {knotU, knotV});
    hold on
    nrbctrlplot(srf12)
    
    
    knotU = srf12.knots{1};
    knotV = srf12.knots{2};
    coefs = srf12.coefs;
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    indexPatch=8;
    
    tempCoord=srf1.coefs;
    
    % P=tempCoord(:,3,3);
    % hold on
    % plot(P(1),P(2),'*k')
    
    %r = 1;
    % %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    %
    
    coefs(1:3,1,1) = tempCoord(1:3,1,3);
    coefs(1:3,1,2) = tempCoord(1:3,2,3);
    coefs(1:3,1,3) = tempCoord(1:3,3,3);
    
    
    
    coefs(1:3,2,1) = [tempCoord(1,1,3);(tempCoord(2,1,3)+L/2)/2; 0];
    coefs(1:3,2,2) = [tempCoord(1,2,3);(tempCoord(2,2,3)+L/2)/2; 0];
    coefs(1:3,2,3) = [tempCoord(1,3,3);(tempCoord(2,3,3)+L/2)/2; 0];
    
    coefs(1:3,3,1) = [tempCoord(1,1,3); L/2; 0];
    coefs(1:3,3,2) = [tempCoord(1,2,3); L/2; 0];
    coefs(1:3,3,3) = [tempCoord(1,3,3); L/2; 0];
    %
    %
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    %
    srf8 = nrbmak(coefs, {knotU, knotV});
    %
    nrbctrlplot(srf8)
    hold on
    %
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
%     knotU = srf8.knots{1};
%     knotV = srf8.knots{2};
%     coefs = srf8.coefs;
    
    tempSrf8=srf8;
 
    tempSrf8 = nrbpermute(tempSrf8,[2,1]);
    %tempSrf8.coefs=tempSrf8.coefs(:,:,end:-1:1);

    knotU = tempSrf8.knots{1};
    knotV = tempSrf8.knots{2};
    coefs = tempSrf8.coefs;

    GIFTmesh{indexPatch} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    indexPatch=2;
    
    tempCoord=srf1.coefs;
    
    P=tempCoord(:,3,1);
    hold on
    plot(P(1),P(2),'*k')
    
    %r = 1;
    % %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) =[tempCoord(1,3,1); -L/2; 0];
    coefs(1:3,1,2) =[tempCoord(1,2,1); -L/2; 0];
    coefs(1:3,1,3) =[tempCoord(1,1,1); -L/2; 0];
    
    coefs(1:3,2,1) = [tempCoord(1,3,1);(tempCoord(2,3,1)-L/2)/2; 0];
    coefs(1:3,2,2) = [tempCoord(1,2,1);(tempCoord(2,2,1)-L/2)/2; 0];
    coefs(1:3,2,3) = [tempCoord(1,1,1);(tempCoord(2,1,1)-L/2)/2; 0];
    
    coefs(1:3,3,1) = tempCoord(1:3,3,1);
    coefs(1:3,3,2) = tempCoord(1:3,2,1);
    coefs(1:3,3,3) = tempCoord(1:3,1,1);
    % %
    % %
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    %
    srf2 = nrbmak(coefs, {knotU, knotV});
    %
    nrbctrlplot(srf2)
    hold on
    %
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
%     knotU = srf2.knots{1};
%     knotV = srf2.knots{2};
%     coefs = srf2.coefs;
    
    tempSrf2=srf2;
 
    tempSrf2 = nrbpermute(tempSrf2,[2,1]);
    tempSrf2.coefs=tempSrf2.coefs(:,end:-1:1,:);

    knotU = tempSrf2.knots{1};
    knotV = tempSrf2.knots{2};
    coefs = tempSrf2.coefs;

    GIFTmesh{indexPatch} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type, 'torus2')
    
    
    GIFTmesh = cell(4,1);
    
    %=== patch 1==================
    patchIndex=1;
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 2;
    
    %a = sqrt(2)/2;
    a = 1/sqrt(2);
    innerRad=0.2;
    outerRad=1;
    
    coefs(1:3,1,1) = [innerRad;0;0];
    coefs(1:3,1,2) = [innerRad*a;innerRad*a;0];
    coefs(1:3,1,3) = [0;innerRad;0];
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = a;
    coefs(4,1,3) = 1;
    
    
    coefs(1:3,2,1) = [outerRad;0;0];
    coefs(1:3,2,2) = [outerRad*a;outerRad*a;0];
    coefs(1:3,2,3) = [0;outerRad;0];
    
    coefs(4,2,1) = 1;
    coefs(4,2,2) = a;
    coefs(4,2,3) = 1;
    
    knotU = [0 0 1 1];
    knotV = [0 0 0 1 1 1];
    surf1 = nrbmak(coefs, {knotU, knotV});
    
    knotU = surf1.knots{1};
    knotV = surf1.knots{2};
    coefs = surf1.coefs;
    
    %     nrbctrlplot(surf1)
    %     hold on
    %     axis equal
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs , p, q, numberElementsU, numberElementsV );
    
    %%%================ patch 2 ===================
    patchIndex=2;
    rotAngle = -pi/2;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotz(rotAngle);
    T2 = vectrans(-pnt);
    
    surf2  = nrbtform(surf1 ,T*rotMatrix*T2);
    
    nrbctrlplot(surf2)
    hold on
    axis equal
    %     pause
    
    knotU = surf2.knots{1};
    knotV = surf2.knots{2};
    coefs = surf2.coefs;
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV,  coefs, p, q, numberElementsU, numberElementsV );
    %%%================ patch 3 ===================
    patchIndex=3;
    rotAngle = -pi/2;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotz(rotAngle);
    T2 = vectrans(-pnt);
    
    surf3  = nrbtform(surf2 ,T*rotMatrix*T2);
    
    %     nrbctrlplot(surf3)
    %     hold on
    %     axis equal
    
    knotU = surf3.knots{1};
    knotV = surf3.knots{2};
    coefs = surf3.coefs;
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV,  coefs, p, q, numberElementsU, numberElementsV );
    
    %%%================ patch 4 ===================
    patchIndex=4;
    rotAngle = -pi/2;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotz(rotAngle);
    T2 = vectrans(-pnt);
    
    surf4  = nrbtform(surf3 ,T*rotMatrix*T2);
    
    %     nrbctrlplot(surf4)
    %     hold on
    %     axis equal
    %
    knotU = surf4.knots{1};
    knotV = surf4.knots{2};
    coefs = surf4.coefs;
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV,  coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%================================
    
elseif strcmp(object_type, 'torus')
    
    GIFTmesh = cell(4,1);
    
    %=== patch 1==================
    patchIndex=1;
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 1;
    
    controlPts = zeros(4,3,2);
    
    innerRad=0.2;
    outerRad=1;
    controlPts(1:4,1,1) =[innerRad,0,0,1];
    controlPts(1:4,2,1) =[innerRad,innerRad,0,1/sqrt(2)];
    controlPts(1:4,3,1) =[0,innerRad,0,1];
    
    controlPts(1:4,1,2) =[outerRad,0,0,1];
    controlPts(1:4,2,2) =[outerRad,outerRad,0,1/sqrt(2)];
    controlPts(1:4,3,2) =[0,outerRad,0,1];
    
    controlPts(1:3,2,1)=controlPts(1:3,2,1)*controlPts(4,2,1);
    controlPts(1:3,2,2)=controlPts(1:3,2,2)*controlPts(4,2,2);
    
    knotU = [0 0 0 1 1 1];
    knotV = [0 0 1 1];
    surf1 = nrbmak(controlPts, {knotU, knotV});
    nrbctrlplot(surf1)
    hold on
    axis equal
    % pause
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, controlPts, p, q, numberElementsU, numberElementsV );
    
    %%%================ patch 2 ===================
    patchIndex=2;
    rotAngle = -pi/2;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotz(rotAngle);
    T2 = vectrans(-pnt);
    
    surf2  = nrbtform(surf1 ,T*rotMatrix*T2);
    
    nrbctrlplot(surf2)
    hold on
    axis equal
    
    
    knotU = surf2.knots{1};
    knotV = surf2.knots{2};
    coefs = surf2.coefs;
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV,  coefs, p, q, numberElementsU, numberElementsV );
    
    %%%================ patch 3 ===================
    patchIndex=3;
    rotAngle = -pi/2;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotz(rotAngle);
    T2 = vectrans(-pnt);
    
    surf3  = nrbtform(surf2 ,T*rotMatrix*T2);
    
    nrbctrlplot(surf3)
    hold on
    axis equal
    
    knotU = surf3.knots{1};
    knotV = surf3.knots{2};
    coefs = surf3.coefs;
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV,  coefs, p, q, numberElementsU, numberElementsV );
    
    %%%================ patch 4 ===================
    patchIndex=4;
    rotAngle = -pi/2;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotz(rotAngle);
    T2 = vectrans(-pnt);
    
    surf4  = nrbtform(surf3 ,T*rotMatrix*T2);
    
    nrbctrlplot(surf4)
    hold on
    axis equal
    
    knotU = surf4.knots{1};
    knotV = surf4.knots{2};
    coefs = surf4.coefs;
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV,  coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%==========================================================================================
elseif strcmp(object_type, '5PatchCircle')
    
    GIFTmesh = cell(numPatches,1);
    %=== patch 1==================
    patchIndex=1;
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    controlPts = zeros(4,3,3);
    alpha=pi/6;
    rad=0.5;
    sq2 = cos(pi/4) ;
    beta = pi/2 - alpha ;
    %sqb = cos(beta/2) ;
    
    rad0 = rad*cos(beta) ;
    rad1 = rad0 + (rad - rad*cos(beta))/2 ;
    
    P1 = [ rad           0               0   1 ;    % pt1
        rad             rad             0   sq2 ;
        0               rad             0   1 ;
        rad1            0               0   1 ;     % pt4
        (rad+rad0/2)/2  (rad+rad0/2)/2  0   sq2 ;
        0               rad1            0   1 ;
        rad0            0               0   1 ;     % pt7
        rad0/2          rad0/2          0   1 ;
        0               rad0            0   1 ] ;
    
    controlPts(1:4,1,1) = P1(1,:)';
    controlPts(1:4,2,1) = P1(2,:)';
    controlPts(1:4,3,1) = P1(3,:)';
    
    controlPts(1:4,1,2) = P1(4,:)';
    controlPts(1:4,2,2) = P1(5,:)';
    controlPts(1:4,3,2) = P1(6,:)';
    
    controlPts(1:4,1,3) = P1(7,:)';
    controlPts(1:4,2,3) = P1(8,:)';
    controlPts(1:4,3,3) = P1(9,:)';
    
    controlPts(1:3,1,1)=controlPts(1:3,1,1)*controlPts(4,1,1);
    controlPts(1:3,2,1) = controlPts(1:3,2,1)*controlPts(4,2,1);
    controlPts(1:3,3,1) = controlPts(1:3,3,1)*controlPts(4,3,1);
    
    controlPts(1:3,1,2) =controlPts(1:3,1,2)*controlPts(4,1,2);
    controlPts(1:3,2,2) =controlPts(1:3,2,2)*controlPts(4,2,2);
    controlPts(1:3,3,2) =controlPts(1:3,3,2)*controlPts(4,3,2);
    
    
    controlPts(1:3,1,3) = controlPts(1:3,1,3)*controlPts(4,1,3);
    controlPts(1:3,2,3) = controlPts(1:3,2,3)*controlPts(4,2,3);
    controlPts(1:3,3,3) = controlPts(1:3,3,3)*controlPts(4,3,3);
    
    knotU = [0 0 0 1 1 1];
    knotV = [0 0 0 1 1 1];
    %     surf = nrbmak(controlPts, {knotU, knotV});
    
    %     nrbctrlplot(surf)
    %     hold on
    %     axis equal
    %     pause
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, controlPts, p, q, numberElementsU, numberElementsV );
    %============================= patch 2 ==========================
    patchIndex=2;
    P2 = P1;
    rotang = pi/2 ; % rotation angle
    ROT = [cos(rotang),sin(rotang),0;-sin(rotang),cos(rotang),0;0,0,1];
    P2(:,1:3) = P2(:,1:3)*ROT ;
    
    controlPts(1:4,1,1) = P2(3,:)';
    controlPts(1:4,2,1) = P2(2,:)';
    controlPts(1:4,3,1) = P2(1,:)';
    
    controlPts(1:4,1,2) = P2(6,:)';
    controlPts(1:4,2,2) = P2(5,:)';
    controlPts(1:4,3,2) = P2(4,:)';
    
    controlPts(1:4,1,3) = P2(9,:)';
    controlPts(1:4,2,3) = P2(8,:)';
    controlPts(1:4,3,3) = P2(7,:)';
    
    controlPts(1:3,1,1)=controlPts(1:3,1,1)*controlPts(4,1,1);
    controlPts(1:3,2,1) = controlPts(1:3,2,1)*controlPts(4,2,1);
    controlPts(1:3,3,1) = controlPts(1:3,3,1)*controlPts(4,3,1);
    
    controlPts(1:3,1,2) =controlPts(1:3,1,2)*controlPts(4,1,2);
    controlPts(1:3,2,2) =controlPts(1:3,2,2)*controlPts(4,2,2);
    controlPts(1:3,3,2) =controlPts(1:3,3,2)*controlPts(4,3,2);
    
    
    controlPts(1:3,1,3) = controlPts(1:3,1,3)*controlPts(4,1,3);
    controlPts(1:3,2,3) = controlPts(1:3,2,3)*controlPts(4,2,3);
    controlPts(1:3,3,3) = controlPts(1:3,3,3)*controlPts(4,3,3);
    
    knotU = [0 0 0 1 1 1];
    knotV = [0 0 0 1 1 1];
    surf = nrbmak(controlPts, {knotU, knotV});
    
    %   nrbctrlplot(surf)
    %     hold on
    %     axis equal
    %     pause
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, controlPts, p, q, numberElementsU, numberElementsV );
    %============================= patch 3 ==========================
    patchIndex=3;
    P3 = P1;
    rotang = pi ; % rotation angle
    ROT = [cos(rotang),sin(rotang),0;-sin(rotang),cos(rotang),0;0,0,1];
    P3(:,1:3) = P3(:,1:3)*ROT ;
    
    controlPts(1:4,1,1) = P3(3,:)';
    controlPts(1:4,2,1) = P3(2,:)';
    controlPts(1:4,3,1) = P3(1,:)';
    
    controlPts(1:4,1,2) = P3(6,:)';
    controlPts(1:4,2,2) = P3(5,:)';
    controlPts(1:4,3,2) = P3(4,:)';
    
    controlPts(1:4,1,3) = P3(9,:)';
    controlPts(1:4,2,3) = P3(8,:)';
    controlPts(1:4,3,3) = P3(7,:)';
    
    controlPts(1:3,1,1)=controlPts(1:3,1,1)*controlPts(4,1,1);
    controlPts(1:3,2,1) = controlPts(1:3,2,1)*controlPts(4,2,1);
    controlPts(1:3,3,1) = controlPts(1:3,3,1)*controlPts(4,3,1);
    
    controlPts(1:3,1,2) =controlPts(1:3,1,2)*controlPts(4,1,2);
    controlPts(1:3,2,2) =controlPts(1:3,2,2)*controlPts(4,2,2);
    controlPts(1:3,3,2) =controlPts(1:3,3,2)*controlPts(4,3,2);
    
    controlPts(1:3,1,3) = controlPts(1:3,1,3)*controlPts(4,1,3);
    controlPts(1:3,2,3) = controlPts(1:3,2,3)*controlPts(4,2,3);
    controlPts(1:3,3,3) = controlPts(1:3,3,3)*controlPts(4,3,3);
    
    knotU = [0 0 0 1 1 1];
    knotV = [0 0 0 1 1 1];
    surf = nrbmak(controlPts, {knotU, knotV});
    
    %nrbctrlplot(surf)
    %hold on
    % axis equal
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, controlPts, p, q, numberElementsU, numberElementsV );
    %============== patch 4 ====================
    patchIndex=4;
    P4 = P1;
    rotang = 3*pi/2 ; % rotation angle
    ROT = [cos(rotang),sin(rotang),0;-sin(rotang),cos(rotang),0;0,0,1];
    P4(:,1:3) = P4(:,1:3)*ROT ;
    
    
    controlPts(1:4,1,1) = P4(1,:)';
    controlPts(1:4,2,1) = P4(2,:)';
    controlPts(1:4,3,1) = P4(3,:)';
    
    controlPts(1:4,1,2) = P4(4,:)';
    controlPts(1:4,2,2) = P4(5,:)';
    controlPts(1:4,3,2) = P4(6,:)';
    
    controlPts(1:4,1,3) = P4(7,:)';
    controlPts(1:4,2,3) = P4(8,:)';
    controlPts(1:4,3,3) = P4(9,:)';
    
    controlPts(1:3,1,1)=controlPts(1:3,1,1)*controlPts(4,1,1);
    controlPts(1:3,2,1) = controlPts(1:3,2,1)*controlPts(4,2,1);
    controlPts(1:3,3,1) = controlPts(1:3,3,1)*controlPts(4,3,1);
    
    controlPts(1:3,1,2) =controlPts(1:3,1,2)*controlPts(4,1,2);
    controlPts(1:3,2,2) =controlPts(1:3,2,2)*controlPts(4,2,2);
    controlPts(1:3,3,2) =controlPts(1:3,3,2)*controlPts(4,3,2);
    
    
    controlPts(1:3,1,3) = controlPts(1:3,1,3)*controlPts(4,1,3);
    controlPts(1:3,2,3) = controlPts(1:3,2,3)*controlPts(4,2,3);
    controlPts(1:3,3,3) = controlPts(1:3,3,3)*controlPts(4,3,3);
    
    knotU = [0 0 0 1 1 1];
    knotV = [0 0 0 1 1 1];
    surf = nrbmak(controlPts, {knotU, knotV});
    
    %nrbctrlplot(surf)
    %hold on
    %axis equal
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, controlPts, p, q, numberElementsU, numberElementsV );
    %============ patch 5 ==================
    
    patchIndex=5;
    P5 = [rad0            0               0   1 ; %
        rad0/2          rad0/2          0   1 ;
        0               rad0            0   1 ;
        rad0/2          -rad0/2         0   1 ; %
        0               0               0   1 ;
        -rad0/2         rad0/2          0   1 ;
        0               -rad0           0   1 ; %
        -rad0/2         -rad0/2         0   1 ;
        -rad0           0               0   1 ];
    
    
    controlPts(1:4,1,1) = P5(7,:)';
    controlPts(1:4,2,1) = P5(4,:)';
    controlPts(1:4,3,1) = P5(1,:)';
    
    controlPts(1:4,1,2) = P5(8,:)';
    controlPts(1:4,2,2) = P5(5,:)';
    controlPts(1:4,3,2) = P5(2,:)';
    
    controlPts(1:4,1,3) = P5(9,:)';
    controlPts(1:4,2,3) = P5(6,:)';
    controlPts(1:4,3,3) = P5(3,:)';
    
    
    controlPts(1:3,1,1)=controlPts(1:3,1,1)*controlPts(4,1,1);
    controlPts(1:3,2,1) = controlPts(1:3,2,1)*controlPts(4,2,1);
    controlPts(1:3,3,1) = controlPts(1:3,3,1)*controlPts(4,3,1);
    
    controlPts(1:3,1,2) =controlPts(1:3,1,2)*controlPts(4,1,2);
    controlPts(1:3,2,2) =controlPts(1:3,2,2)*controlPts(4,2,2);
    controlPts(1:3,3,2) =controlPts(1:3,3,2)*controlPts(4,3,2);
    
    
    controlPts(1:3,1,3) = controlPts(1:3,1,3)*controlPts(4,1,3);
    controlPts(1:3,2,3) = controlPts(1:3,2,3)*controlPts(4,2,3);
    controlPts(1:3,3,3) = controlPts(1:3,3,3)*controlPts(4,3,3);
    
    knotU = [0 0 0 1 1 1];
    knotV = [0 0 0 1 1 1];
    surf = nrbmak(controlPts, {knotU, knotV});
    
    %nrbctrlplot(surf)
    %hold on
    %axis equal
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, controlPts, p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type, '5patchSquare')
    
    GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    
    % ======================= patch 1 =========================
    patchIndex=1;
    
    coefs(1:3,1,1) = [0; 0; 0]/4;
    coefs(1:3,1,2) = [1; 1; 0]/4;
    coefs(1:3,2,1) = [4; 0; 0]/4;
    coefs(1:3,2,2) = [3; 1; 0]/4;
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    
    % ======================= patch 2 =========================
    patchIndex=2;
    
    coefs(1:3,1,1) = [4; 0; 0]/4;
    coefs(1:3,1,2) = [3; 1; 0]/4;
    coefs(1:3,2,1) = [4; 4; 0]/4;
    coefs(1:3,2,2) = [3; 3; 0]/4;
    
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    % ======================= patch 3 =========================
    patchIndex=3;
    
    coefs(1:3,1,1) = [0; 4; 0]/4;
    coefs(1:3,1,2) = [1; 3; 0]/4;
    coefs(1:3,2,1) = [4; 4; 0]/4;
    coefs(1:3,2,2) = [3; 3; 0]/4;
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    % ======================= patch 4 =========================
    patchIndex=4;
    
    coefs(1:3,1,1) = [0; 0; 0]/4;
    coefs(1:3,1,2) = [1; 1; 0]/4;
    coefs(1:3,2,1) = [0; 4; 0]/4;
    coefs(1:3,2,2) = [1; 3; 0]/4;
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    % ======================= patch 5 =========================
    
    patchIndex=5;
    
    coefs(1:3,1,1) = [1; 1; 0]/4;
    coefs(1:3,1,2) = [1; 3; 0]/4;
    coefs(1:3,2,1) = [3; 1; 0]/4;
    coefs(1:3,2,2) = [3; 3; 0]/4;
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    
    GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    %=========================================================================================================
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
elseif strcmp(object_type, 'NURBSsphere')
    
    % patch 1  bottom ==================================
    patchIndex=1;
    numberElementsU = 1;
    numberElementsV = 1;
    p = 4;
    q = 4;
    
    row0=[4*(1-sqrt(3)), 4*(1-sqrt(3)),4*(1-sqrt(3)),4*(3-sqrt(3));...
        -sqrt(2),sqrt(2)*(sqrt(3)-4),sqrt(2)*(sqrt(3)-4),sqrt(2)*(3*sqrt(3)-2);...
        0, 4*(1-2*sqrt(3))/3,4*(1-2*sqrt(3))/3,4*(5-sqrt(3))/3;...
        sqrt(2),sqrt(2)*(sqrt(3)-4),sqrt(2)*(sqrt(3)-4),sqrt(2)*(3*sqrt(3)-2);...
        4*(sqrt(3)-1),4*(1-sqrt(3)),4*(1-sqrt(3)),4*(3-sqrt(3))];
    
    row1=[sqrt(2)*(sqrt(3)-4),-sqrt(2),sqrt(2)*(sqrt(3)-4),sqrt(2)*(3*sqrt(3)-2);...
        (2-3*sqrt(3))/2,(2-3*sqrt(3))/2,-(sqrt(3)+6)/2,(sqrt(3)+6)/2;...
        0,sqrt(2)*(2*sqrt(3)-7)/3,-5*sqrt(6)/3,sqrt(2)*(sqrt(3)+6)/3;...
        (3*sqrt(3)-2)/2, (2-3*sqrt(3))/2,-(sqrt(3)+6)/2,(sqrt(3)+6)/2;...
        sqrt(2)*(4-sqrt(3)),-sqrt(2),sqrt(2)*(sqrt(3)-4),sqrt(2)*(3*sqrt(3)-2)];
    
    row2=[4*(1-2*sqrt(3))/3, 0,4*(1-2*sqrt(3))/3,4*(5-sqrt(3))/3;...
        sqrt(2)*(2*sqrt(3)-7)/3, 0, -5*sqrt(6)/3,sqrt(2)*(sqrt(3)+6)/3;...
        0, 0, 4*(sqrt(3)-5)/3,4*(5*sqrt(3)-1)/9;...
        sqrt(2)*(7-2*sqrt(3))/3,0,-5*sqrt(6)/3,sqrt(2)*(sqrt(3)+6)/3;...
        4*(2*sqrt(3)-1)/3,0,4*(1-2*sqrt(3))/3,4*(5-sqrt(3))/3];
    
    row3=[sqrt(2)*(sqrt(3)-4),sqrt(2),sqrt(2)*(sqrt(3)-4),sqrt(2)*(3*sqrt(3)-2);...
        (2-3*sqrt(3))/2,(3*sqrt(3)-2)/2,-(sqrt(3)+6)/2,(sqrt(3)+6)/2;...
        0, sqrt(2)*(7-2*sqrt(3))/3,-5*sqrt(6)/3,sqrt(2)*(sqrt(3)+6)/3;...
        (3*sqrt(3)-2)/2,(3*sqrt(3)-2)/2,-(sqrt(3)+6)/2,(sqrt(3)+6)/2;...
        sqrt(2)*(4-sqrt(3)),sqrt(2),sqrt(2)*(sqrt(3)-4),sqrt(2)*(3*sqrt(3)-2)];
    
    row4=[4*(1-sqrt(3)),4*(sqrt(3)-1),4*(1-sqrt(3)),4*(3-sqrt(3));...
        -sqrt(2),sqrt(2)*(4-sqrt(3)),sqrt(2)*(sqrt(3)-4),sqrt(2)*(3*sqrt(3)-2);...
        0, 4*(2*sqrt(3)-1)/3,4*(1-2*sqrt(3))/3,4*(5-sqrt(3))/3;...
        sqrt(2),sqrt(2)*(4-sqrt(3)),sqrt(2)*(sqrt(3)-4),sqrt(2)*(3*sqrt(3)-2);...
        4*(sqrt(3)-1), 4*(sqrt(3)-1),  4*(1-sqrt(3)), 4*(3-sqrt(3))];
    
    controlPts(1:4,1,1) = row0(1,:);
    controlPts(1:4,2,1) =  row0(2,:);
    controlPts(1:4,3,1) =  row0(3,:);
    controlPts(1:4,4,1) =  row0(4,:);
    controlPts(1:4,5,1) =  row0(5,:);
    
    controlPts(1:4,1,2) = row1(1,:);
    controlPts(1:4,2,2) =row1(2,:);
    controlPts(1:4,3,2) =row1(3,:);
    controlPts(1:4,4,2) =row1(4,:);
    controlPts(1:4,5,2) =row1(5,:);
    
    controlPts(1:4,1,3) =row2(1,:);
    controlPts(1:4,2,3) =row2(2,:);
    controlPts(1:4,3,3) =row2(3,:);
    controlPts(1:4,4,3) =row2(4,:);
    controlPts(1:4,5,3) =row2(5,:);
    
    controlPts(1:4,1,4) =row3(1,:);
    controlPts(1:4,2,4) =row3(2,:);
    controlPts(1:4,3,4) =row3(3,:);
    controlPts(1:4,4,4) =row3(4,:);
    controlPts(1:4,5,4) =row3(5,:);
    
    controlPts(1:4,1,5) =row4(1,:);
    controlPts(1:4,2,5) =row4(2,:);
    controlPts(1:4,3,5) =row4(3,:);
    controlPts(1:4,4,5) =row4(4,:);
    controlPts(1:4,5,5) =row4(5,:);
    
    knotU = [0 0 0 0 0 1 1 1 1 1];
    knotV = [0 0 0 0 0 1 1 1 1 1];
    surf1 = nrbmak(controlPts, {knotU, knotV});
    
    knotU = surf1.knots{1};
    knotV = surf1.knots{2};
    coefs = surf1.coefs;
    %GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    GIFTmesh{patchIndex} = genGIFTmesh_thinShell( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %     nrbctrlplot(surf1)
    %     pause
    %     hold on
    
    
    % patch 4 front ===========================
    
    rotAngle = pi/2;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotx(rotAngle);
    T2 = vectrans(-pnt);
    
    surf2  = nrbtform(surf1 ,T*rotMatrix*T2);
    
    rotAngle = pi/2;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotx(rotAngle);
    T2 = vectrans(-pnt);
    
    surf3  = nrbtform(surf2 ,T*rotMatrix*T2);
    
    patchIndex=4;
    rotAngle = pi/2;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotx(rotAngle);
    T2 = vectrans(-pnt);
    
    surf4  = nrbtform(surf3 ,T*rotMatrix*T2);
    
    knotU = surf4.knots{1};
    knotV = surf4.knots{2};
    coefs = surf4.coefs;
    GIFTmesh{patchIndex} = genGIFTmesh_thinShell( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %     nrbkntplot(surf4 )
    %     hold on
    %     axis equal
    
    % patch 5 left ===========================
    patchIndex=5;
    rotAngle = pi/2;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecroty(rotAngle);
    T2 = vectrans(-pnt);
    
    surf5  = nrbtform(surf1 ,T*rotMatrix*T2);
    
    surf5 = nrbpermute(surf5,[2,1]);
    
    knotU = surf5.knots{1};
    knotV = surf5.knots{2};
    coefs = surf5.coefs;
    GIFTmesh{patchIndex} = genGIFTmesh_thinShell( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %     nrbkntplot(surf5 )
    %     hold on
    %     axis equal
    
    % patch 6 rigth ===========================
    patchIndex=6;
    rotAngle = pi/2;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotz(rotAngle);
    T2 = vectrans(-pnt);
    surf6  = nrbtform(surf4 ,T*rotMatrix*T2);
    
    
    knotU = surf6.knots{1};
    knotV = surf6.knots{2};
    coefs = surf6.coefs;
    
    GIFTmesh{patchIndex} = genGIFTmesh_thinShell( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %     nrbkntplot(surf6 )
    %     hold on
    %     axis equal
    % patch 3 up ==================================
    
    rotAngle = pi/2;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotz(rotAngle);
    T2 = vectrans(-pnt);
    surf2  = nrbtform(surf6 ,T*rotMatrix*T2);
    patchIndex=3;
    rotAngle = pi/2;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotx(rotAngle);
    T2 = vectrans(-pnt);
    
    surf3  = nrbtform(surf2 ,T*rotMatrix*T2);
    
    knotU = surf3.knots{1};
    knotV = surf3.knots{2};
    coefs = surf3.coefs;
    tempCoefs=coefs;
    tempCoefs(:,1,1)=coefs(:,5,1);
    tempCoefs(:,2,1)=coefs(:,4,1);
    tempCoefs(:,4,1)=coefs(:,2,1);
    tempCoefs(:,5,1)=coefs(:,1,1);
    
    tempCoefs(:,1,2)=coefs(:,5,2);
    tempCoefs(:,2,2)=coefs(:,4,2);
    tempCoefs(:,4,2)=coefs(:,2,2);
    tempCoefs(:,5,2)=coefs(:,1,2);
    
    tempCoefs(:,1,3)=coefs(:,5,3);
    tempCoefs(:,2,3)=coefs(:,4,3);
    tempCoefs(:,4,3)=coefs(:,2,3);
    tempCoefs(:,5,3)=coefs(:,1,3);
    
    tempCoefs(:,1,4)=coefs(:,5,4);
    tempCoefs(:,2,4)=coefs(:,4,4);
    tempCoefs(:,4,4)=coefs(:,2,4);
    tempCoefs(:,5,4)=coefs(:,1,4);
    
    tempCoefs(:,1,5)=coefs(:,5,5);
    tempCoefs(:,2,5)=coefs(:,4,5);
    tempCoefs(:,4,5)=coefs(:,2,5);
    tempCoefs(:,5,5)=coefs(:,1,5);
    
    coefs = tempCoefs;
    surf3.coefs= tempCoefs;
    GIFTmesh{patchIndex} = genGIFTmesh_thinShell( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %     nrbkntplot(surf3 )
    %     hold on
    %     axis equal
    %============patch 2 back =========
    patchIndex=2;
    
    rotAngle = -pi/2;
    pnt = [0,0,0];
    T  = vectrans(pnt);
    rotMatrix = vecrotz(rotAngle);
    T2 = vectrans(-pnt);
    
    surf2  = nrbtform(surf5 ,T*rotMatrix*T2);
    
    knotU = surf2.knots{1};
    knotV = surf2.knots{2};
    coefs = surf2.coefs;
    %GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    GIFTmesh{patchIndex} = genGIFTmesh_thinShell( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %     nrbkntplot(surf2 )
    %     hold on
    %     axis equal
    
    
    tempGIFTmesh=GIFTmesh;
    GIFTmesh{1}=tempGIFTmesh{4};
    GIFTmesh{2}=tempGIFTmesh{6};
    GIFTmesh{3}=tempGIFTmesh{2};
    GIFTmesh{4}=tempGIFTmesh{5};
    GIFTmesh{5}=tempGIFTmesh{1};
    GIFTmesh{6}=tempGIFTmesh{3};
    
elseif strcmp(object_type, 'pinchedCylinder')
    
    R      = 300;
    L      = 300;
    z = L;
    
    GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 1;
    
    for patchIndex = 1:numPatches
        
        %initialize geometry on coarsest mesh
        controlPts = zeros(4,3,2);
        
        % first zKnot
        
        controlPts(1:3,1,1) = [R;0;0];
        controlPts(1:3,2,1) = [R;R;0];
        controlPts(1:3,3,1) = [0;R;0];
        
        % third zKnot
        z = L;
        controlPts(1:3,1,2) = [R;0;z];
        controlPts(1:3,2,2) = [R;R;z];
        controlPts(1:3,3,2) = [0;R;z];
        
        controlPts(4,:,:)   = 1;
        
        fac                 = 1/sqrt(2);
        
        controlPts(4,2,1) = fac;
        controlPts(4,2,2) = fac;
        
        % homogenous coordinates (x*w,y*w,z*w)
        
        controlPts(1:3,2,1) = controlPts(1:3,2,1)*fac;
        controlPts(1:3,2,2) = controlPts(1:3,2,2)*fac;
        
        knotU = [0 0 0 1 1 1];
        knotV = [0 0 1 1];
        surf = nrbmak(controlPts, {knotU, knotV})
        
        nrbctrlplot(surf)
        
        GIFTmesh{patchIndex} = genGIFTmesh_thinShell( knotU, knotV, controlPts, p, q, numberElementsU, numberElementsV );
        
    end
    
elseif strcmp(object_type, 'pinchedCylinder2patch')
    
    GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 1;
    controlPts = zeros(4,3,2);
    
    %data for the first path
    controlPts(:,:,1) = 1e2*[3.000000000000000   2.560660171779821   1.810660171779821
        0   1.060660171779821   1.810660171779821
        0                   0                   0
        0.010000000000000   0.008535533905933   0.008535533905933];
    
    controlPts(:,:,2) = 1e2*[3.000000000000000   2.560660171779821   1.810660171779821
        0   1.060660171779821   1.810660171779821
        3.000000000000000   2.560660171779821   2.560660171779821
        0.010000000000000   0.008535533905933   0.008535533905933 ];
    
    knotU = [0 0 0 1 1 1];
    knotV = [0 0 1 1];
    
    surf1 = nrbmak(controlPts, {knotU, knotV})
    nrbctrlplot(surf1)
    
    GIFTmesh{1} = genGIFTmesh_thinShell( knotU, knotV, controlPts, p, q, numberElementsU, numberElementsV );
    
    %data for the second patch
    
    %initialize geometry on coarsest mesh
    controlPts = zeros(4,3,2);
    
    controlPts(:,:,1) = 1e2*[   1.810660171779821   1.060660171779821                   0
        1.810660171779821   2.560660171779821   3.000000000000000
        0                   0                   0
        0.008535533905933   0.008535533905933   0.010000000000000];
    
    controlPts(:,:,2) = 1e2*[1.810660171779821   1.060660171779821                   0
        1.810660171779821   2.560660171779821   3.000000000000000
        2.560660171779821   2.560660171779821   3.000000000000000
        0.008535533905933   0.008535533905933   0.010000000000000];
    
    knotU = [0 0 0 1 1 1];
    knotV = [0 0 1 1];
    surf2 = nrbmak(controlPts, {knotU, knotV})
    
    hold on
    nrbctrlplot(surf2)
    
    GIFTmesh{2} = genGIFTmesh_thinShell( knotU, knotV, controlPts, p, q, numberElementsU, numberElementsV );
    
    
elseif strcmp(object_type, 'pipeShell')
    
    
elseif strcmp(object_type, 'pinchedHalfCylinder2patch')
    
    radius=300;
    height=600;
    
    R=radius;
    z=height;
    
    numPatches = 2;
    GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    
    p = 2;
    q = 1;
    
    patchIndex=1;
    coefs = zeros(4,3,2);
    
    coefs(1:3,1,1) = [R;0;0];  %in
    coefs(1:3,2,1) = [R;R;0];  %in
    coefs(1:3,3,1) = [0;R;0];  %in
    
    coefs(1:3,1,2) = [R;0;z];  %in
    coefs(1:3,2,2) = [R;R;z];  %in
    coefs(1:3,3,2) = [0;R;z];  %in
    
    fac                 = 1/sqrt(2);
    coefs(4,:,:)   = 1;
    coefs(4,2,1) = fac;  %in
    coefs(4,2,2) = fac;  %in
    
    coefs(1:3,2,1) = coefs(1:3,2,1)*fac;  %in
    coefs(1:3,2,2) = coefs(1:3,2,2)*fac;  %in
    
    knotU = [0 0 0 1 1 1];
    knotV = [0 0 1 1];
    surf1 = nrbmak(coefs, {knotU, knotV});
    
    nrbctrlplot(surf1)
    
    GIFTmesh{patchIndex} = genGIFTmesh_thinShell( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    patchIndex=2;
    coefs = zeros(4,3,2);
    
    coefs(1:3,1,1) = [0;R;0];
    coefs(1:3,2,1) = [-R;R;0];
    coefs(1:3,3,1) = [-R;0;0];
    
    coefs(1:3,1,2) = [0;R;z];
    coefs(1:3,2,2) = [-R;R;z];
    coefs(1:3,3,2) = [-R;0;z];
    
    fac                 = 1/sqrt(2);
    coefs(4,:,:)   = 1;
    coefs(4,2,1) = fac;  %in
    coefs(4,2,2) = fac;  %in
    
    coefs(1:3,2,1) = coefs(1:3,2,1)*fac;  %in
    coefs(1:3,2,2) = coefs(1:3,2,2)*fac;  %in
    
    knotU = [0 0 0 1 1 1];
    knotV = [0 0 1 1];
    surf2 = nrbmak(coefs, {knotU, knotV});
    hold on
    
    nrbctrlplot(surf2)
    
    GIFTmesh{patchIndex} = genGIFTmesh_thinShell( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    
elseif strcmp(object_type, 'plate_hole')
    
    %initialize two-patch geometry
    GIFTmesh = cell(2,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %data for the first patch%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 1;
    
    %initialize geometry on coarsest mesh
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    rad = 1;
    side_fac = 1;
    coefs(1:3,1,1) = rad.*[-1;0;0];
    coefs(1:3,2,1) = rad.*[-0.853553390593274; 0.353553390593274; 0];
    coefs(1:3,3,1) = rad.*[-0.603553390593274; 0.603553390593274; 0];
    coefs(1:3,1,2) = side_fac.*[-4;0;0];
    coefs(1:3,2,2) = side_fac.*[-4;2;0];
    coefs(1:3,3,2) = side_fac.*[-4;4;0];
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 0.853553390593274;
    coefs(4,3,1) = 0.853553390593274;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    coefs(4,3,2) = 1;
    
    GIFTmesh{1} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % data for second patch%
    %%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 1;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    coefs(1:3,1,1) = rad.*[-0.603553390593274;0.603553390593274;0];
    coefs(1:3,2,1) = rad.*[-0.353553390593274;0.853553390593274;0];
    coefs(1:3,3,1) = rad.*[0;1;0];
    coefs(1:3,1,2) = side_fac.*[-4;4;0];
    coefs(1:3,2,2) = side_fac.*[-2;4;0];
    coefs(1:3,3,2) = side_fac.*[0;4;0];
    coefs(4,1,1) = 0.853553390593274;
    coefs(4,2,1) = 0.853553390593274;
    coefs(4,3,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    coefs(4,3,2) = 1;
    GIFTmesh{2} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type, 'Lshape_trap')
    %L-shaped domain with exact solution discretized with 2 trapezoidal patches
    numPatches = 2;
    GIFTmesh = cell(numPatches,1);
    
    %%%%%%%%%%%%%
    %first patch%
    %%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    coefs(1:3,1,1) = [-1/sqrt(2); 1/sqrt(2); 0];
    coefs(1:3,1,2) = [0; sqrt(2); 0];
    coefs(1:3,2,1) = [0; 0; 0];
    coefs(1:3,2,2) = [sqrt(2); 0; 0];
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    GIFTmesh{1} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%
    %second patch%
    %%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    coefs(1:3,1,1) = [0; 0; 0];
    coefs(1:3,1,2) = [sqrt(2); 0; 0];
    coefs(1:3,2,1) = [-1/sqrt(2); -1/sqrt(2); 0];
    coefs(1:3,2,2) = [0; -sqrt(2); 0];
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    GIFTmesh{2} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    
elseif strcmp(object_type, 'plate_hole2')
    
    %initialize two-patch geometry
    GIFTmesh = cell(2,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %data for the first patch%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %initialize geometry on coarsest mesh
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    rad = 1;
    side_fac = 1;
    coefs(1:3,1,1) = rad.*[-1;0;0];
    coefs(1:3,2,1) = rad.*[-0.603553390593274; 0.603553390593274; 0];
    coefs(1:3,1,2) = side_fac.*[-4;0;0];
    coefs(1:3,2,2) = side_fac.*[-4;4;0];
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    
    GIFTmesh{1} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % data for second patch%
    %%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    coefs(1:3,1,1) = rad.*[-0.603553390593274;0.603553390593274;0];
    coefs(1:3,2,1) = rad.*[0;1;0];
    coefs(1:3,1,2) = side_fac.*[-4;4;0];
    coefs(1:3,2,2) = side_fac.*[0;4;0];
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    GIFTmesh{2} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type, 'bilinear_wedge')
    
    %initialize two-patch geometry
    GIFTmesh = cell(2,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %data for the first patch%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %initialize geometry on coarsest mesh
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    rad = 1;
    side_fac = 1;
    coefs(1:3,1,1) = [-3;-1;0];
    coefs(1:3,2,1) = [0; 0; 0];
    coefs(1:3,1,2) = [-2;2;0];
    coefs(1:3,2,2) = [0;1;0];
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    
    GIFTmesh{1} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % data for second patch%
    %%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    coefs(1:3,1,1) = [0;0;0];
    coefs(1:3,2,1) = [3;1/2;0];
    coefs(1:3,1,2) = [0;1;0];
    coefs(1:3,2,2) = [5/2;3;0];
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    GIFTmesh{2} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type, 'bilinear_rectangle')
    
    %initialize two-patch geometry
    GIFTmesh = cell(2,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %data for the first patch%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %initialize geometry on coarsest mesh
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    rad = 1;
    side_fac = 1;
    coefs(1:3,1,1) = [-1;0;0];
    coefs(1:3,2,1) = [0;0;0];
    coefs(1:3,1,2) = [-1;1;0];
    coefs(1:3,2,2) = [0;1;0];
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    
    GIFTmesh{1} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % data for second patch%
    %%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    coefs(1:3,1,1) = [0;0;0];
    coefs(1:3,2,1) = [1;0;0];
    coefs(1:3,1,2) = [0;1;0];
    %coefs(1:3,2,2) = [1;1.2;0];
    coefs(1:3,2,2) = [1;1;0];
    % coefs(1:3,2,2) = [1;1.5;0];
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    GIFTmesh{2} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type, 'bilinear_rectangleTest')
    
    %initialize two-patch geometry
    GIFTmesh = cell(2,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %data for the first patch%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %initialize geometry on coarsest mesh
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    rad = 1;
    side_fac = 1;
    coefs(1:3,1,1) = [-1;0;0];
    coefs(1:3,2,1) = [0; 0; 0];
    coefs(1:3,1,2) = [-1;1;0];
    coefs(1:3,2,2) = [0;1;0];
    
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    
    GIFTmesh{1} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % data for second patch%
    %%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    coefs(1:3,1,1) = [0;0;0];
    coefs(1:3,2,1) = [1;0;0];
    coefs(1:3,1,2) = [0;1;0];
    coefs(1:3,2,2) = [1;1;0];
    
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    
    solid1 = nrbmak(coefs,{knotU knotV});
    solid1 = nrbpermute(solid1,[2,1]);
    knotU=solid1.knots{1};
    knotV=solid1.knots{2};
    coefs=solid1.coefs;
    
    GIFTmesh{2} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    %-----
elseif strcmp(object_type, 'threePatch')
    
    %initialize two-patch geometry
    GIFTmesh = cell(3,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %data for the first patch%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %initialize geometry on coarsest mesh
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    rad = 1;
    side_fac = 1;
    %     coefs(1:3,1,1) = [0;0;0];
    %     coefs(1:3,2,1) = [2;0; 0];
    %     coefs(1:3,1,2) = [1;1;0];
    %     coefs(1:3,2,2) = [2;1;0];
    
    coefs(1:3,1,1) = [0;0;0]/2;
    coefs(1:3,2,1) = [2;0; 0]/2;
    coefs(1:3,1,2) = [1;1;0]/2;
    coefs(1:3,2,2) = [2;1;0]/2;
    
    
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    
    GIFTmesh{1} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % data for second patch%
    %%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    %     coefs(1:3,1,1) = [1;1;0];
    %     coefs(1:3,2,1) = [2;1;0];
    %     coefs(1:3,1,2) = [1;2;0];
    %     coefs(1:3,2,2) = [2;2;0];
    
    coefs(1:3,1,1) = [1;1;0]/2;
    coefs(1:3,2,1) = [2;1;0]/2;
    coefs(1:3,1,2) = [1;2;0]/2;
    coefs(1:3,2,2) = [2;2;0]/2;
    
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    
    %     solid1 = nrbmak(coefs,{knotU knotV});
    %     solid1 = nrbpermute(solid1,[2,1]);
    %     knotU=solid1.knots{1};
    %     knotV=solid1.knots{2};
    %     coefs=solid1.coefs;
    
    GIFTmesh{2} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % data for 3rd patch%
    %%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    %     coefs(1:3,1,1) = [0;0;0];
    %     coefs(1:3,2,1) = [1;1;0];
    %     coefs(1:3,1,2) = [0;2;0];
    %     coefs(1:3,2,2) = [1;2;0];
    
    coefs(1:3,1,1) = [0;0;0]/2;
    coefs(1:3,2,1) = [1;1;0]/2;
    coefs(1:3,1,2) = [0;2;0]/2;
    coefs(1:3,2,2) = [1;2;0]/2;
    
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    
    %     solid1 = nrbmak(coefs,{knotU knotV});
    %     solid1 = nrbpermute(solid1,[2,1]);
    %     knotU=solid1.knots{1};
    %     knotV=solid1.knots{2};
    %     coefs=solid1.coefs;
    
    GIFTmesh{3} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type, 'bilinear_rectangleUpDown')
    
    %initialize two-patch geometry
    GIFTmesh = cell(2,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %data for the first patch%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %initialize geometry on coarsest mesh
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    rad = 1;
    side_fac = 1;
    coefs(1:3,1,1) = [-1;0;0];
    coefs(1:3,2,1) = [0; 0; 0];
    coefs(1:3,1,2) = [-1;1;0];
    coefs(1:3,2,2) = [0;1;0];
    %
    %     coefs(1:3,1,2) = [-1;0;0];
    %     coefs(1:3,2,2) = [0; 0; 0];
    %     coefs(1:3,1,1) = [-1;1;0];
    %     coefs(1:3,2,1) = [0;1;0];
    
    
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    
    GIFTmesh{1} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % data for second patch%
    %%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    %     coefs(1:3,1,1) = [-1;1;0];
    %     coefs(1:3,2,1) = [0;1;0];
    %     coefs(1:3,1,2) = [-1;2;0];
    %     coefs(1:3,2,2) = [0;2;0];
    
    coefs(1:3,1,2) = [-1;1;0];
    coefs(1:3,2,2) = [0;1;0];
    coefs(1:3,1,1) = [-1;2;0];
    coefs(1:3,2,1) = [0;2;0];
    
    
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    
    solid1 = nrbmak(coefs,{knotU knotV});
    solid1 = nrbpermute(solid1,[2,1]);
    knotU=solid1.knots{1};
    knotV=solid1.knots{2};
    coefs=solid1.coefs;
    
    GIFTmesh{2} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    
elseif strcmp(object_type, 'nonas_rectangle')
    
    %initialize two-patch geometry
    GIFTmesh = cell(2,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %data for the first patch%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %initialize geometry on coarsest mesh
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    coefs(1:3,1,1) = [-1;0;0];
    coefs(1:3,2,1) = [0; 0; 0];
    coefs(1:3,1,2) = [-1;1;0];
    coefs(1:3,2,2) = [0;1;0];
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    
    GIFTmesh{1} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % data for second patch%
    %%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 2;
    p = 2;
    q = 1;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0.5, 1, 1];
    
    %initialize geometry on coarsest mesh
    coefs(1:3,1,1) = [0;0;0];
    coefs(1:3,2,1) = [0.25;0;0];
    coefs(1:3,3,1) = [1;0;0];
    coefs(1:3,1,2) = [0;.5;0];
    coefs(1:3,2,2) = [0.5;0.5;0];
    coefs(1:3,3,2) = [1;0.5;0];
    coefs(1:3,1,3) = [0;1;0];
    coefs(1:3,2,3) = [0.25;1;0];
    coefs(1:3,3,3) = [1;1;0];
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 1;
    coefs(4,3,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    coefs(4,3,2) = 1;
    coefs(4,1,3) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,3) = 1;
    GIFTmesh{2} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type, 'nonas_rectangle_curved')
    
    %initialize two-patch geometry
    GIFTmesh = cell(2,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %data for the first patch%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %initialize geometry on coarsest mesh
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    coefs(1:3,1,1) = [-1;0;0];
    coefs(1:3,2,1) = [0; 0; 0];
    coefs(1:3,1,2) = [-1;1;0];
    coefs(1:3,2,2) = [0;1;0];
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    
    GIFTmesh{1} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % data for second patch%
    %%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 2;
    numberElementsV = 2;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 0.5, 1, 1, 1];
    knotV = [0, 0, 0, 0.5, 1, 1, 1];
    
    %initialize geometry on coarsest mesh
    coefs(1:3,1,1) = [0;0;0];
    coefs(1:3,2,1) = [0.25;0;0];
    coefs(1:3,3,1) = [0.75;0;0];
    coefs(1:3,4,1) = [1;0;0];
    coefs(1:3,1,2) = [0;.25;0];
    coefs(1:3,2,2) = [0.5;0.25;0];
    coefs(1:3,3,2) = [0.9;0.25;0];
    coefs(1:3,4,2) = [1;0.25;0];
    coefs(1:3,1,3) = [0;.75;0];
    coefs(1:3,2,3) = [0.5;0.75;0];
    coefs(1:3,3,3) = [0.9;0.75;0];
    coefs(1:3,4,3) = [1;0.75;0];
    coefs(1:3,1,4) = [0;1;0];
    coefs(1:3,2,4) = [0.25;1;0];
    coefs(1:3,3,4) = [0.75;1;0];
    coefs(1:3,4,4) = [1;1;0];
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 1;
    coefs(4,3,1) = 1;
    coefs(4,4,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    coefs(4,3,2) = 1;
    coefs(4,4,2) = 1;
    coefs(4,1,3) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,3) = 1;
    coefs(4,4,3) = 1;
    coefs(4,1,4) = 1;
    coefs(4,2,4) = 1;
    coefs(4,3,4) = 1;
    coefs(4,4,4) = 1;
    GIFTmesh{2} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type, 'bilinear_vShape')
    
    %initialize two-patch geometry
    GIFTmesh = cell(2,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %data for the first patch%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %initialize geometry on coarsest mesh
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    rad = 1;
    side_fac = 1;
    coefs(1:3,1,1) = [-1;1;0];
    coefs(1:3,2,1) = [0; 0; 0];
    coefs(1:3,1,2) = [-1;2;0];
    coefs(1:3,2,2) = [0;1;0];
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    
    GIFTmesh{1} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % data for second patch%
    %%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    coefs(1:3,1,1) = [0;0;0];
    coefs(1:3,2,1) = [1;1;0];
    coefs(1:3,1,2) = [0;1;0];
    coefs(1:3,2,2) = [1;2;0];
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    GIFTmesh{2} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type, 'bilinear_nonGeneric')
    
    %initialize two-patch geometry
    GIFTmesh = cell(2,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %data for the first patch%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %initialize geometry on coarsest mesh
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    rad = 1;
    side_fac = 1;
    coefs(1:3,1,1) = [-3;-1;0];
    coefs(1:3,2,1) = [0;0;0];
    coefs(1:3,1,2) = [-2;2;0];
    coefs(1:3,2,2) = [0;1;0];
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    
    GIFTmesh{1} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % data for second patch%
    %%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    coefs(1:3,1,1) = [0;0;0];
    coefs(1:3,2,1) = [3;1/2;0];
    coefs(1:3,1,2) = [0;1;0];
    coefs(1:3,2,2) = [5/2;3;0];
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    GIFTmesh{2} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type, 'plate_holeC1')
    numberElementsU = 2;
    numberElementsV = 1;
    p = 2;
    q = 2;
    knotU = [0, 0, 0, 0.5, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    weights = [1,(1+1/sqrt(2))/2,(1+1/sqrt(2))/2,1, 1,1,1,1, 1,1,1,1]';
    controlPoints = [-1, 0; -1, sqrt(2)-1; 1-sqrt(2), 1; 0, 1; -2.5, 0; -2.5, 0.75; -0.75, 2.5; 0, 2.5; -4, 0; -4, 4; -4, 4; 0, 4];
    %the number of control points in the u and v directions
    lenU = length(knotU)-3;  %number of basis functions in the  u direction
    lenV = length(knotV)-3;  %number of basis functions in the v direction
    coordinates = [controlPoints, weights];
    
    %convert control points/weights to coefs format in NURBS toolbox
    coefs = zeros(4,lenU,lenV);
    for i=1:lenU
        for j=1:lenV
            index = (j-1)*lenU+i;
            xcoord = coordinates(index,1);
            ycoord = coordinates(index,2);
            wght = coordinates(index,3);
            coefs(1,i,j) = xcoord*wght;
            coefs(2,i,j) = ycoord*wght;
            coefs(4,i,j) = wght;
        end
    end
    GIFTmesh{1} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
elseif strcmp(object_type, 'Lshape_ex')
    %L-shaped domain with exact solution discretized with 3 patches
    numPatches = 3;
    GIFTmesh = cell(numPatches,1);
    
    %%%%%%%%%%%%%
    %first patch%
    %%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    coefs(1:3,1,1) = [0; 0; 0];
    coefs(1:3,1,2) = [-1/sqrt(2); 1/sqrt(2); 0];
    coefs(1:3,2,1) = [1/sqrt(2); 1/sqrt(2); 0];
    coefs(1:3,2,2) = [0; sqrt(2); 0];
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    GIFTmesh{1} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%
    %second patch%
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    coefs(1:3,1,1) = [0; 0; 0];
    coefs(1:3,1,2) = [1/sqrt(2); 1/sqrt(2); 0];
    coefs(1:3,2,1) = [1/sqrt(2); -1/sqrt(2); 0];
    coefs(1:3,2,2) = [sqrt(2); 0; 0];
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    GIFTmesh{2} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%
    %third patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    coefs(1:3,1,1) = [0; 0; 0];
    coefs(1:3,1,2) = [1/sqrt(2); -1/sqrt(2); 0];
    coefs(1:3,2,1) = [-1/sqrt(2); -1/sqrt(2); 0];
    coefs(1:3,2,2) = [0; -sqrt(2); 0];
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    GIFTmesh{3} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
elseif strcmp(object_type, 'EdgeCrackMultiPatch')
    
    GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %divide the patches along the x direction
    xVertices = linspace(0,L,3);
    yVertices = linspace(0,W,3);
    
    patchCounter = 0;
    
    for patchIndexY = 1:(numPatches-2)
        for patchIndexX = 1:(numPatches-2)
            %set the dimensions of the patch
            patchCounter = patchCounter + 1;
            patchMinX = xVertices(patchIndexX);
            patchMaxX = xVertices(patchIndexX+1);
            patchMinY = yVertices(patchIndexY);
            patchMaxY = yVertices(patchIndexY+1);
            
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
            
            GIFTmesh{patchCounter} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
            
            
        end
    end
    
    tempGIFTmesh = GIFTmesh{4};
    GIFTmesh{4} = GIFTmesh{3};
    GIFTmesh{3} = tempGIFTmesh;
elseif strcmp(object_type, 'spanner')
    p=2;
    q=2;
    numberElementsU = 3;
    numberElementsV = 1;
    load('spanner.mat')
    knotU = surf.knots{1};
    knotV = [0,0,0,1,1,1];
    scale_fac = 516/12;
    surf.coefs(1:3,:,:) = scale_fac*surf.coefs(1:3,:,:);
    GIFTmesh{1} = genGIFTmesh( knotU, knotV, surf.coefs(:,:,1:3), p, q, numberElementsU, numberElementsV );
    GIFTmesh{2} = genGIFTmesh( knotU, knotV, surf.coefs(:,:,3:5), p, q, numberElementsU, numberElementsV );
    GIFTmesh{3} = genGIFTmesh( knotU, knotV, surf.coefs(:,:,5:7), p, q, numberElementsU, numberElementsV );
    GIFTmesh{4} = genGIFTmesh( knotU, knotV, surf.coefs(:,:,7:9), p, q, numberElementsU, numberElementsV );
    GIFTmesh{5} = genGIFTmesh( knotU, knotV, surf.coefs(:,:,9:11), p, q, numberElementsU, numberElementsV );
    GIFTmesh{6} = genGIFTmesh( knotU, knotV, surf.coefs(:,:,11:13), p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type, 'circle')
    p=2;
    q=2;
    numberElementsU = 1;
    numberElementsV = 1;
    knotU = [0,0,0,1,1,1];
    knotV = [0,0,0,1,1,1];
    coefs = [-sqrt(2)/4, -sqrt(2)/2, -sqrt(2)/4, 0, 0, 0, sqrt(2)/4, sqrt(2)/2, sqrt(2)/4; sqrt(2)/4, 0, -sqrt(2)/4, sqrt(2)/2, 0, -sqrt(2)/2, sqrt(2)/4, 0, -sqrt(2)/4; zeros(1,9);...
        1, sqrt(2)/2, 1, sqrt(2)/2, 1, sqrt(2)/2, 1, sqrt(2)/2, 1];
    coefs(1,:) = coefs(1,:).*coefs(4,:);
    coefs(2,:) = coefs(2,:).*coefs(4,:);
    coefs = reshape(coefs,[4,3,3]);
    GIFTmesh{1} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type,'LShaped_bracket')
    
    %Multipatch bracket with 4 patches
    %numPatches = 10;
    numPatches = 18;
    GIFTmesh = cell(numPatches,1);
    
    %%%%%%%%%%%%%
    %first patch%
    %%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    a = sqrt(2)/2;
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [a; a; 0];
    coefs(1:3,1,2) = [0; 1; 0];
    coefs(1:3,1,3) = [-a; a; 0];
    coefs(1:3,2,1) = [2; 2; 0];
    coefs(1:3,2,2) = [0; 4; 0];
    coefs(1:3,2,3) = [-2; 2; 0];
    coefs(1:3,3,1) = [4; 4; 0];
    coefs(1:3,3,2) = [0; 4; 0];
    coefs(1:3,3,3) = [-4; 4; 0];
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = a;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    GIFTmesh{1} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    
    %%%%%%%%%%%%%%
    %second patch%
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [-a; a; 0];
    coefs(1:3,1,2) = [-1; 0; 0];
    coefs(1:3,1,3) = [-a; -a; 0];
    coefs(1:3,2,1) = [-2; 2; 0];
    coefs(1:3,2,2) = [-4; 0; 0];
    coefs(1:3,2,3) = [-2; -2; 0];
    coefs(1:3,3,1) = [-4; 4; 0];
    coefs(1:3,3,2) = [-4; 0; 0];
    coefs(1:3,3,3) = [-4; -4; 0];
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = a;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    GIFTmesh{2} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    %%%%%%%%%%%%%%
    %third patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [-a; -a; 0];
    coefs(1:3,1,2) = [0; -1; 0];
    coefs(1:3,1,3) = [a; -a; 0];
    coefs(1:3,2,1) = [-2; -2; 0];
    coefs(1:3,2,2) = [0; -4; 0];
    coefs(1:3,2,3) = [2; -2; 0];
    coefs(1:3,3,1) = [-4; -4; 0];
    coefs(1:3,3,2) = [0; -4; 0];
    coefs(1:3,3,3) = [4; -4; 0];
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = a;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    GIFTmesh{3} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%
    %Fourth patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [a;-a; 0];
    coefs(1:3,1,2) = [1; 0; 0];
    coefs(1:3,1,3) = [a; a; 0];
    coefs(1:3,2,1) = [2; -2; 0];
    coefs(1:3,2,2) = [4; 0; 0];
    coefs(1:3,2,3) = [2; 2; 0];
    coefs(1:3,3,1) = [4; -4; 0];
    coefs(1:3,3,2) = [4; 0; 0];
    coefs(1:3,3,3) = [4; 4; 0];
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = a;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    GIFTmesh{4} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    
    %%%%%%%%%%%%%%
    %Fifth patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [a; a; 0];
    coefs(1:3,1,2) = [2; 2; 0];
    coefs(1:3,1,3) = [4; 4; 0];
    coefs(1:3,2,1) = [1; 0; 0];
    coefs(1:3,2,2) = [4; 0; 0];
    coefs(1:3,2,3) = [4; 0; 0];
    coefs(1:3,3,1) = [a; -a; 0];
    coefs(1:3,3,2) = [2; -2; 0];
    coefs(1:3,3,3) = [4; -4; 0];
    
    coefs(1,:,:,:) = coefs(1,:,:,:)-8;
    coefs(1,2,1) = (1/a-8)*a; %correct the translation because the coefs are using projective coordinates
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = a;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    surf = nrbmak(coefs,{knotU, knotV});
    surf = nrbpermute(surf,[2,1]);
    
    GIFTmesh{5} = genGIFTmesh( surf.knots{1}, surf.knots{2}, surf.coefs, p, q, numberElementsU, numberElementsV );
    
    
    %%%%%%%%%%%%%%
    %Sixth patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [a; -a; 0];
    coefs(1:3,1,2) = [2; -2; 0];
    coefs(1:3,1,3) = [4; -4; 0];
    coefs(1:3,2,1) = [0; -1; 0];
    coefs(1:3,2,2) = [0; -4; 0];
    coefs(1:3,2,3) = [0; -4; 0];
    coefs(1:3,3,1) = [-a; -a; 0];
    coefs(1:3,3,2) = [-2; -2; 0];
    coefs(1:3,3,3) = [-4; -4; 0];
    
    
    coefs(1,:,:,:) = coefs(1,:,:,:)-8;
    coefs(1,2,1) = (((1/a-8)*a)-1); %correct the translation because the coefs are using projective coordinates
    %coefs(2,2,1) = a;
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = a;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    
    %GIFTmesh{6} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    surf = nrbmak(coefs,{knotU, knotV});
    surf = nrbpermute(surf,[2,1]);
    
    GIFTmesh{6} = genGIFTmesh( surf.knots{1}, surf.knots{2}, surf.coefs, p, q, numberElementsU, numberElementsV );
    
    
    %%%%%%%%%%%%%%
    %Seventh patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    
    coefs(1:3,1,1) = [-a; -a; 0];
    coefs(1:3,1,2) = [-2; -2; 0];
    coefs(1:3,1,3) = [-4; -4; 0];
    coefs(1:3,2,1) = [-1; 0; 0];
    coefs(1:3,2,2) = [-4; 0; 0];
    coefs(1:3,2,3) = [-4; 0; 0];
    coefs(1:3,3,1) = [-a; a; 0];
    coefs(1:3,3,2) = [-2; 2; 0];
    coefs(1:3,3,3) = [-4; 4; 0];
    
    
    coefs(1,:,:,:) = coefs(1,:,:,:)-8;
    coefs(1,2,1) = (1/a-8); %correct the translation because the coefs are using projective coordinates
    %coefs(2,2,1) = a;
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = a;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    
    %GIFTmesh{7} = genGIFTmesh(knotU,  knotV, coefs, p, q, numberElementsU, numberElementsV );
    surf = nrbmak(coefs,{knotU, knotV});
    surf = nrbpermute(surf,[2,1]);
    
    GIFTmesh{7} = genGIFTmesh( surf.knots{1}, surf.knots{2}, surf.coefs, p, q, numberElementsU, numberElementsV );
    
    
    
    %%%%%%%%%%%%%%
    %Eigth patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    
    coefs(1:3,1,1) = [-a; a; 0];
    coefs(1:3,1,2) = [-2; 2; 0];
    coefs(1:3,1,3) = [-4; 4; 0];
    coefs(1:3,2,1) = [0; 1; 0];
    coefs(1:3,2,2) = [0; 4; 0];
    coefs(1:3,2,3) = [0; 4; 0];
    coefs(1:3,3,1) = [a; a; 0];
    coefs(1:3,3,2) = [2; 2; 0];
    coefs(1:3,3,3) = [4; 4; 0];
    
    
    coefs(1,:,:,:) = coefs(1,:,:,:)-8;
    coefs(1,2,1) = (((1/a-8)*a)-1); %correct the translation because the coefs are using projective coordinates
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = a;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    %GIFTmesh{8} = genGIFTmesh(knotU,knotV,coefs, p, q, numberElementsU, numberElementsV );
    surf = nrbmak(coefs,{knotU, knotV});
    surf = nrbpermute(surf,[2,1]);
    
    GIFTmesh{8} = genGIFTmesh( surf.knots{1}, surf.knots{2}, surf.coefs, p, q, numberElementsU, numberElementsV );
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %9th patch%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 1;
    
    %initialize geometry on coarsest mesh
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    coefs(1:3,1,1) = [-4;0;0];
    coefs(1:3,2,1) = [-0.853553390593274*4; 0.353553390593274*4; 0];
    coefs(1:3,3,1) = [-0.603553390593274*4;0.603553390593274*4;0];
    coefs(1:3,1,2) = [-12;0;0];
    coefs(1:3,2,2) = [-12;8;0];
    coefs(1:3,3,2) = [-12;12;0];
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 0.853553390593274;
    coefs(4,3,1) = 0.853553390593274;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    coefs(4,3,2) = 1;
    
    trans = vectrans([-12.0 -8.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    
    GIFTmesh{9} =  genGIFTmesh( knotU, knotV, solid2.coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % 10th patch%
    %%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 1;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    
    coefs(1:3,1,1) = [-0.603553390593274*4;0.603553390593274*4;0];
    coefs(1:3,2,1) = [-0.353553390593274*4;0.853553390593274*4;0];
    coefs(1:3,3,1) = [0;4;0];
    coefs(1:3,1,2) = [-12;12;0];
    coefs(1:3,2,2) = [-8;12;0];
    coefs(1:3,3,2) = [0;12;0];
    coefs(4,1,1) = 0.853553390593274;
    coefs(4,2,1) = 0.853553390593274;
    coefs(4,3,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    coefs(4,3,2) = 1;
    
    trans = vectrans([-12.0 -8.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    
    GIFTmesh{10} =  genGIFTmesh( knotU, knotV, solid2.coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%
    %11th patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    a = sqrt(2)/2;
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [a; a; 0];
    coefs(1:3,1,2) = [0; 1; 0];
    coefs(1:3,1,3) = [-a; a; 0];
    coefs(1:3,2,1) = [2; 2; 0];
    coefs(1:3,2,2) = [0; 4; 0];
    coefs(1:3,2,3) = [-2; 2; 0];
    coefs(1:3,3,1) = [4; 4; 0];
    coefs(1:3,3,2) = [0; 4; 0];
    coefs(1:3,3,3) = [-4; 4; 0];
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = a;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    %coefs(1,:,:,:) = coefs(1,:,:,:)-20;
    %coefs(2,:,:,:) = coefs(2,:,:,:)-12;
    %coefs(1,1,2) = -20; %correct the translation because the coefs are using projective coordinates
    %coefs(1,2,2) = (1/a-12);
    
    trans = vectrans([-20.0 -12.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    
    GIFTmesh{11} = genGIFTmesh( knotU, knotV, solid2.coefs, p, q, numberElementsU, numberElementsV );
    
    %patch 12
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [-a; a; 0];
    coefs(1:3,1,2) = [-1; 0; 0];
    coefs(1:3,1,3) = [-a; -a; 0];
    coefs(1:3,2,1) = [-2; 2; 0];
    coefs(1:3,2,2) = [-4; 0; 0];
    coefs(1:3,2,3) = [-2; -2; 0];
    coefs(1:3,3,1) = [-4; 4; 0];
    coefs(1:3,3,2) = [-4; 0; 0];
    coefs(1:3,3,3) = [-4; -4; 0];
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = a;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    trans = vectrans([-20.0 -12.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    
    GIFTmesh{12} = genGIFTmesh( knotU, knotV, solid2.coefs, p, q, numberElementsU, numberElementsV  );
    
    %%%%%%%%%%%%%%
    %13th patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [-a; -a; 0];
    coefs(1:3,1,2) = [0; -1; 0];
    coefs(1:3,1,3) = [a; -a; 0];
    coefs(1:3,2,1) = [-2; -2; 0];
    coefs(1:3,2,2) = [0; -4; 0];
    coefs(1:3,2,3) = [2; -2; 0];
    coefs(1:3,3,1) = [-4; -4; 0];
    coefs(1:3,3,2) = [0; -4; 0];
    coefs(1:3,3,3) = [4; -4; 0];
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = a;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    trans = vectrans([-20.0 -12.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    
    GIFTmesh{13} = genGIFTmesh( knotU, knotV, solid2.coefs, p, q, numberElementsU, numberElementsV  );
    
    %%%%%%%%%%%%%%
    %14th patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [a;-a; 0];
    coefs(1:3,1,2) = [1; 0; 0];
    coefs(1:3,1,3) = [a; a; 0];
    coefs(1:3,2,1) = [2; -2; 0];
    coefs(1:3,2,2) = [4; 0; 0];
    coefs(1:3,2,3) = [2; 2; 0];
    coefs(1:3,3,1) = [4; -4; 0];
    coefs(1:3,3,2) = [4; 0; 0];
    coefs(1:3,3,3) = [4; 4; 0];
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = a;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    trans = vectrans([-20.0 -12.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    
    GIFTmesh{14} = genGIFTmesh( knotU, knotV, solid2.coefs, p, q, numberElementsU, numberElementsV  );
    
    
    
    %%%%%%%%%%%%%%
    %15th patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [-a; a; 0];
    coefs(1:3,1,2) = [-2; 2; 0];
    coefs(1:3,1,3) = [-4; 4; 0];
    coefs(1:3,2,1) = [0; 1; 0];
    coefs(1:3,2,2) = [0; 4; 0];
    coefs(1:3,2,3) = [0; 4; 0];
    coefs(1:3,3,1) = [a; a; 0];
    coefs(1:3,3,2) = [2; 2; 0];
    coefs(1:3,3,3) = [4; 4; 0];
    
    coefs(1,:,:,:) = coefs(1,:,:,:)-8;
    coefs(1,2,1) = (((1/a-8)*a)-1); %correct the translation because the coefs are using projective coordinates
    
    %coefs(1,2,1) = (1/a-8)*a; %correct the translation because the coefs are using projective coordinates
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = a;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    trans = vectrans([-12.0 -20.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    
    %GIFTmesh{15} = genGIFTmesh( knotU, knotV, solid2.coefs, p, q, numberElementsU, numberElementsV );
    surf = nrbmak(solid2.coefs,{knotU, knotV});
    surf = nrbpermute(surf,[2,1]);
    
    GIFTmesh{15} = genGIFTmesh( surf.knots{1}, surf.knots{2}, surf.coefs, p, q, numberElementsU, numberElementsV );
    
    
    %%%%%%%%%%%%%%
    %16th patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [-a; -a; 0];
    coefs(1:3,1,2) = [-2; -2; 0];
    coefs(1:3,1,3) = [-4; -4; 0];
    coefs(1:3,2,1) = [-1; 0; 0];
    coefs(1:3,2,2) = [-4; 0; 0];
    coefs(1:3,2,3) = [-4; 0; 0];
    coefs(1:3,3,1) = [-a; a; 0];
    coefs(1:3,3,2) = [-2; 2; 0];
    coefs(1:3,3,3) = [-4; 4; 0];
    
    coefs(1,:,:,:) = coefs(1,:,:,:)-8;
    %coefs(1,2,1) = (((1/a-8)*a)-1); %correct the translation because the coefs are using projective coordinates
    coefs(1,2,1) = (1/a-8); %correct the translation because the coefs are using projective coordinates
    
    %coefs(2,2,1) = a;
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = a;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    
    trans = vectrans([-12.0 -20.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    
    
    surf = nrbmak(solid2.coefs,{knotU, knotV});
    surf = nrbpermute(surf,[2,1]);
    
    GIFTmesh{16} = genGIFTmesh( surf.knots{1}, surf.knots{2}, surf.coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%
    %17th patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [a; -a; 0];
    coefs(1:3,1,2) = [2; -2; 0];
    coefs(1:3,1,3) = [4; -4; 0];
    coefs(1:3,2,1) = [0; -1; 0];
    coefs(1:3,2,2) = [0; -4; 0];
    coefs(1:3,2,3) = [0; -4; 0];
    coefs(1:3,3,1) = [-a; -a; 0];
    coefs(1:3,3,2) = [-2; -2; 0];
    coefs(1:3,3,3) = [-4; -4; 0];
    
    coefs(1,:,:,:) = coefs(1,:,:,:)-8;
    coefs(1,2,1) = (((1/a-8)*a)-1); %correct the translation because the coefs are using projective coordinates
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = a;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    
    trans = vectrans([-12.0 -20.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    
    surf = nrbmak(solid2.coefs,{knotU, knotV});
    surf = nrbpermute(surf,[2,1]);
    
    GIFTmesh{17} = genGIFTmesh( surf.knots{1}, surf.knots{2}, surf.coefs, p, q, numberElementsU, numberElementsV );
    
    
    %GIFTmesh{17} = genGIFTmesh( knotU, knotV, solid2.coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%
    %18th patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [a; a; 0];
    coefs(1:3,1,2) = [2; 2; 0];
    coefs(1:3,1,3) = [4; 4; 0];
    coefs(1:3,2,1) = [1; 0; 0];
    coefs(1:3,2,2) = [4; 0; 0];
    coefs(1:3,2,3) = [4; 0; 0];
    coefs(1:3,3,1) = [a; -a; 0];
    coefs(1:3,3,2) = [2; -2; 0];
    coefs(1:3,3,3) = [4; -4; 0];
    
    coefs(1,:,:,:) = coefs(1,:,:,:)-8;
    coefs(1,2,1) = (1/a-8)*a; %correct the translation because the coefs are using projective coordinates
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = a;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    trans = vectrans([-12.0 -20.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    
    surf = nrbmak(solid2.coefs,{knotU, knotV});
    surf = nrbpermute(surf,[2,1]);
    
    GIFTmesh{18} = genGIFTmesh( surf.knots{1}, surf.knots{2}, surf.coefs, p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type,'LShaped_bracketTest')
    
    %Multipatch bracket with 4 patches
    %numPatches = 10;
    numPatches = 2;
    GIFTmesh = cell(numPatches,1);
    
    %%%%%%%%%%%%%
    %first patch%
    %%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    a = sqrt(2)/2;
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    %     coefs(1:3,1,1) = [a; a; 0];
    %     coefs(1:3,1,2) = [0; 1; 0];
    %     coefs(1:3,1,3) = [-a; a; 0];
    %     coefs(1:3,2,1) = [2; 2; 0];
    %     coefs(1:3,2,2) = [0; 4; 0];
    %     coefs(1:3,2,3) = [-2; 2; 0];
    %     coefs(1:3,3,1) = [4; 4; 0];
    %     coefs(1:3,3,2) = [0; 4; 0];
    %     coefs(1:3,3,3) = [-4; 4; 0];
    %
    %     coefs(4,1,1) = 1;
    %     coefs(4,1,2) = a;
    %     coefs(4,1,3) = 1;
    %     coefs(4,2,1) = 1;
    %     coefs(4,2,2) = 1;
    %     coefs(4,2,3) = 1;
    %     coefs(4,3,1) = 1;
    %     coefs(4,3,2) = 1;
    %     coefs(4,3,3) = 1;
    
    coefs(1:3,3,3) = [a; a; 0];
    coefs(1:3,3,2) = [0; 1; 0];
    coefs(1:3,3,1) = [-a; a; 0];
    coefs(1:3,2,3) = [2; 2; 0];
    coefs(1:3,2,2) = [0; 4; 0];
    coefs(1:3,2,1) = [-2; 2; 0];
    coefs(1:3,1,3) = [4; 4; 0];
    coefs(1:3,1,2) = [0; 4; 0];
    coefs(1:3,1,1) = [-4; 4; 0];
    
    coefs(4,3,3) = 1;
    coefs(4,3,2) = a;
    coefs(4,3,1) = 1;
    coefs(4,2,3) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,1,3) = 1;
    coefs(4,1,2) = 1;
    coefs(4,1,1) = 1;
    
    
    GIFTmesh{2} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    
    %%%%%%%%%%%%%%
    %second patch%
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    %     coefs(1:3,1,1) = [-a; a; 0];
    %     coefs(1:3,1,2) = [-1; 0; 0];
    %     coefs(1:3,1,3) = [-a; -a; 0];
    %     coefs(1:3,2,1) = [-2; 2; 0];
    %     coefs(1:3,2,2) = [-4; 0; 0];
    %     coefs(1:3,2,3) = [-2; -2; 0];
    %     coefs(1:3,3,1) = [-4; 4; 0];
    %     coefs(1:3,3,2) = [-4; 0; 0];
    %     coefs(1:3,3,3) = [-4; -4; 0];
    %
    %     coefs(4,1,1) = 1;
    %     coefs(4,1,2) = a;
    %     coefs(4,1,3) = 1;
    %     coefs(4,2,1) = 1;
    %     coefs(4,2,2) = 1;
    %     coefs(4,2,3) = 1;
    %     coefs(4,3,1) = 1;
    %     coefs(4,3,2) = 1;
    %     coefs(4,3,3) = 1;
    
    coefs(1:3,3,3) = [-a; a; 0];
    coefs(1:3,3,2) = [-1; 0; 0];
    coefs(1:3,3,1) = [-a; -a; 0];
    coefs(1:3,2,3) = [-2; 2; 0];
    coefs(1:3,2,2) = [-4; 0; 0];
    coefs(1:3,2,1) = [-2; -2; 0];
    coefs(1:3,1,3) = [-4; 4; 0];
    coefs(1:3,1,2) = [-4; 0; 0];
    coefs(1:3,1,1) = [-4; -4; 0];
    
    coefs(4,3,3) = 1;
    coefs(4,3,2) = a;
    coefs(4,3,1) = 1;
    coefs(4,2,3) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,1,3) = 1;
    coefs(4,1,2) = 1;
    coefs(4,1,1) = 1;
    
    solid1 = nrbmak(coefs,{knotU knotV});
    solid1 = nrbpermute(solid1,[2,1]);
    
    knotU=solid1.knots{1};
    knotV=solid1.knots{2};
    coefs=solid1.coefs;
    
    GIFTmesh{1} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
end
