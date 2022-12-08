function [PHUTelemV,GIFTmeshV,PHUTelemCR,GIFTmeshCR,PHUTelemCL,GIFTmeshCL] = rotateGiftMeshAndPHUTelem2(numberElementsU,numberElementsV)
%rotate patches on 5patches example

%============================== v shape ===================================
pointV =[ -6.8284    2.8284;...
    -5.4142    2.8284;...
    -4.0000         0;...
    -4.0000    1.4142;...
    -1.1716    2.8284;...
    -2.5858    2.8284];

coefsVA(1:3,1,1)=[pointV(1,:),0];
coefsVA(1:3,1,2)=[pointV(2,:),0];
coefsVA(1:3,2,1) =[pointV(3,:),0];
coefsVA(1:3,2,2) =[pointV(4,:),0];

coefsVB(1:3,1,1)=[pointV(3,:),0];
coefsVB(1:3,1,2)=[pointV(4,:),0];
coefsVB(1:3,2,1) =[pointV(5,:),0];
coefsVB(1:3,2,2) =[pointV(6,:),0];

coefsVA(4,1,1) = 1;
coefsVA(4,1,2) = 1;
coefsVA(4,2,1) = 1;
coefsVA(4,2,2) = 1;

coefsVB(4,1,1) = 1;
coefsVB(4,1,2) = 1;
coefsVB(4,2,1) = 1;
coefsVB(4,2,2) = 1;

knotU = [0 0 1 1];
knotV = [0 0 1 1];
p = 1;
q = 1;
GIFTmeshElementsU=1;
GIFTmeshElementsV=1;

GIFTmeshV{1}= genGIFTmesh( knotU, knotV, coefsVA, p, q, GIFTmeshElementsU, GIFTmeshElementsV );
[PHUTelemV{1}, ~] = initPHTmeshGen(p,q,numberElementsU,numberElementsV);

GIFTmeshV{2}= genGIFTmesh( knotU, knotV, coefsVB, p, q, GIFTmeshElementsU, GIFTmeshElementsV );
[PHUTelemV{2}, ~] = initPHTmeshGen(p,q,numberElementsU,numberElementsV);

% figure
% plotPHTMeshMP2(PHUTelemV{1}, GIFTmeshV{1});
% hold on
% plotPHTMeshMP2(PHUTelemV{2}, GIFTmeshV{2});
% axis equal
% pause

%============================ center right ====================================

pointCR =[-0.0000   -2.0000;
    0.0000    2.0000;
    -1.0000   -1.0000;
    -1.0000    1.0000;
    -3.0000   -1.0000;
    -3.0000    1.0000];

coefsCRB(1:3,1,1)=[pointCR(3,:),0];
coefsCRB(1:3,1,2)=[pointCR(4,:),0];
coefsCRB(1:3,2,1) =[pointCR(1,:),0];
coefsCRB(1:3,2,2) =[pointCR(2,:),0];

coefsCRA(1:3,1,1) =[pointCR(5,:),0];
coefsCRA(1:3,1,2) =[pointCR(6,:),0];
coefsCRA(1:3,2,1)=[pointCR(3,:),0];
coefsCRA(1:3,2,2)=[pointCR(4,:),0];

% coefsCRA(1:3,1,1)=[pointCR(3,:),0];
% coefsCRA(1:3,1,2)=[pointCR(4,:),0];
% coefsCRA(1:3,2,1) =[pointCR(1,:),0];
% coefsCRA(1:3,2,2) =[pointCR(2,:),0];
% 
% coefsCRB(1:3,1,1) =[pointCR(5,:),0];
% coefsCRB(1:3,1,2) =[pointCR(6,:),0];
% coefsCRB(1:3,2,1)=[pointCR(3,:),0];
% coefsCRB(1:3,2,2)=[pointCR(4,:),0];


coefsCRA(4,1,1) = 1;
coefsCRA(4,1,2) = 1;
coefsCRA(4,2,1) = 1;
coefsCRA(4,2,2) = 1;

coefsCRB(4,1,1) = 1;
coefsCRB(4,1,2) = 1;
coefsCRB(4,2,1) = 1;
coefsCRB(4,2,2) = 1;

knotU = [0 0 1 1];
knotV = [0 0 1 1];
p = 1;
q = 1;
GIFTmeshElementsU=1;
GIFTmeshElementsV=1;

GIFTmeshCR{1}= genGIFTmesh( knotU, knotV, coefsCRA, p, q, GIFTmeshElementsU, GIFTmeshElementsV );
[PHUTelemCR{1}, ~] = initPHTmeshGen(p,q,numberElementsU,numberElementsV);

GIFTmeshCR{2}= genGIFTmesh( knotU, knotV, coefsCRB, p, q, GIFTmeshElementsU, GIFTmeshElementsV );
[PHUTelemCR{2}, ~] = initPHTmeshGen(p,q,numberElementsU,numberElementsV);

% figure
% plotPHTMeshMP2(PHUTelemCR{1}, GIFTmeshCR{1});
% hold on
% pause
% plotPHTMeshMP2(PHUTelemCR{2}, GIFTmeshCR{2});
% axis equal
% pause
%============== center left ==================

pointCL =[-2    -2;
    -2     2;
    -1    -1;
    -1     1;
    1    -1;
    1     1];

coefsCLA(1:3,1,1)=[pointCL(1,:),0];
coefsCLA(1:3,1,2)=[pointCL(2,:),0];
coefsCLA(1:3,2,1) =[pointCL(3,:),0];
coefsCLA(1:3,2,2) =[pointCL(4,:),0];

coefsCLB(1:3,1,1)=[pointCL(3,:),0];
coefsCLB(1:3,1,2)=[pointCL(4,:),0];
coefsCLB(1:3,2,1) =[pointCL(5,:),0];
coefsCLB(1:3,2,2) =[pointCL(6,:),0];

coefsCLA(4,1,1) = 1;
coefsCLA(4,1,2) = 1;
coefsCLA(4,2,1) = 1;
coefsCLA(4,2,2) = 1;

coefsCLB(4,1,1) = 1;
coefsCLB(4,1,2) = 1;
coefsCLB(4,2,1) = 1;
coefsCLB(4,2,2) = 1;

knotU = [0 0 1 1];
knotV = [0 0 1 1];
p = 1;
q = 1;
GIFTmeshElementsU=1;
GIFTmeshElementsV=1;

GIFTmeshCL{1}= genGIFTmesh( knotU, knotV, coefsCLA, p, q, GIFTmeshElementsU, GIFTmeshElementsV );
[PHUTelemCL{1}, ~] = initPHTmeshGen(p,q,numberElementsU,numberElementsV);

GIFTmeshCL{2}= genGIFTmesh( knotU, knotV, coefsCLB, p, q, GIFTmeshElementsU, GIFTmeshElementsV );
[PHUTelemCL{2}, ~] = initPHTmeshGen(p,q,numberElementsU,numberElementsV);

% figure
% plotPHTMeshMP2(PHUTelemCL{1}, GIFTmeshCL{1});
% hold on
% pause
% plotPHTMeshMP2(PHUTelemCL{2}, GIFTmeshCL{2});
% axis equal

end

