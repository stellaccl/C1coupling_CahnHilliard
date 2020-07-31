function [PHUTelemTBleft,GIFTmeshTBleft,PHUTelemTBright,GIFTmeshTBright] = rotateGiftMeshAndPHUTelemTopBottom(numberElementsU,numberElementsV)
%rotate patches on 5patches example

%============================== top bottom shape ===================================
% coefsB(1:3,1,1) = [2; 0; 0];
% coefsB(1:3,1,2) = [1; 1; 0];
% coefsB(1:3,2,1) = [2; 4; 0];
% coefsB(1:3,2,2) = [1; 3; 0];

% coefsA(1:3,1,1) = [0; 0; 0];
% coefsA(1:3,1,2) = [1; 1; 0];
% coefsA(1:3,2,1) = [0; 4; 0];
% coefsA(1:3,2,2) = [1; 3; 0];

coefsB(1:3,1,1) = [4; 0; 0];
coefsB(1:3,1,2) = [4; 4; 0];
coefsB(1:3,2,1) = [8; 0; 0];
coefsB(1:3,2,2) = [8; 4; 0];

coefsA(1:3,1,1) = [3; 1; 0];
coefsA(1:3,1,2) = [3; 3; 0];
coefsA(1:3,2,1) = [4; 0; 0];
coefsA(1:3,2,2) = [4; 4; 0];

coefsA(4,1,1) = 1;
coefsA(4,1,2) = 1;
coefsA(4,2,1) = 1;
coefsA(4,2,2) = 1;

coefsB(4,1,1) = 1;
coefsB(4,1,2) = 1;
coefsB(4,2,1) = 1;
coefsB(4,2,2) = 1;

knotU = [0 0 1 1];
knotV = [0 0 1 1];
p = 1;
q = 1;
GIFTmeshElementsU=1;
GIFTmeshElementsV=1;

GIFTmeshTBleft{1}= genGIFTmesh( knotU, knotV, coefsA, p, q, GIFTmeshElementsU, GIFTmeshElementsV );
[PHUTelemTBleft{1}, ~] = initPHTmeshGen(p,q,numberElementsU,numberElementsV);

GIFTmeshTBleft{2}= genGIFTmesh( knotU, knotV, coefsB, p, q, GIFTmeshElementsU, GIFTmeshElementsV );
[PHUTelemTBleft{2}, ~] = initPHTmeshGen(p,q,numberElementsU,numberElementsV);


% figure
% plotPHTMeshMP(PHUTelemTBleft, GIFTmeshTBleft);
% pause

coefsA(1:3,1,1) = [8; 0; 0];
coefsA(1:3,1,2) = [8; 4; 0];
coefsA(1:3,2,1) = [9; 1; 0];
coefsA(1:3,2,2) = [9; 3; 0];

coefsB(1:3,1,1) = [4; 0; 0];
coefsB(1:3,1,2) = [4; 4; 0];
coefsB(1:3,2,1) = [8; 0; 0];
coefsB(1:3,2,2) = [8; 4; 0];

coefsA(4,1,1) = 1;
coefsA(4,1,2) = 1;
coefsA(4,2,1) = 1;
coefsA(4,2,2) = 1;

coefsB(4,1,1) = 1;
coefsB(4,1,2) = 1;
coefsB(4,2,1) = 1;
coefsB(4,2,2) = 1;

knotU = [0 0 1 1];
knotV = [0 0 1 1];
p = 1;
q = 1;
GIFTmeshElementsU=1;
GIFTmeshElementsV=1;

GIFTmeshTBright{1}= genGIFTmesh( knotU, knotV, coefsA, p, q, GIFTmeshElementsU, GIFTmeshElementsV );
[PHUTelemTBright{1}, ~] = initPHTmeshGen(p,q,numberElementsU,numberElementsV);

GIFTmeshTBright{2}= genGIFTmesh( knotU, knotV, coefsB, p, q, GIFTmeshElementsU, GIFTmeshElementsV );
[PHUTelemTBright{2}, ~] = initPHTmeshGen(p,q,numberElementsU,numberElementsV);

% figure
% plotPHTMeshMP(PHUTelemTBright, GIFTmeshTBright);
% pause
end

