function [a,b,c] = computeFunctionABC_GIFTmesh_uDirection(PHUTelem,GIFTmesh,patchA,elemA,patchB,elemB,refPoint)

xmin = PHUTelem{patchA}(elemA).vertex(1);
xmax = PHUTelem{patchA}(elemA).vertex(3);
ymin = PHUTelem{patchA}(elemA).vertex(2);
ymax = PHUTelem{patchA}(elemA).vertex(4);

urefA=refPoint.urefA;
vrefA=refPoint.vrefA;
urefB=refPoint.urefB;
vrefB=refPoint.vrefB;

[coord, dxdxi] = paramMap( GIFTmesh{patchA}, urefA, vrefA, xmin, ymin, xmax, ymax);
plot(coord(1),coord(2),'+r')
hold on

dGdx_Patch1=dxdxi(1,:);
dGdy_Patch1=dxdxi(2,:);

xmin = PHUTelem{patchB}(elemB).vertex(1);
xmax = PHUTelem{patchB}(elemB).vertex(3);
ymin = PHUTelem{patchB}(elemB).vertex(2);
ymax = PHUTelem{patchB}(elemB).vertex(4);

[coord, dxdxi] = paramMap( GIFTmesh{patchB}, urefB, vrefB, xmin, ymin, xmax, ymax);
plot(coord(1),coord(2),'.k')
hold on

dGdx_Patch2=dxdxi(1,:);
dGdy_Patch2=dxdxi(2,:);

a=dGdy_Patch1(1)*dGdy_Patch2(2)-dGdy_Patch2(1)*dGdy_Patch1(2);
b=-(dGdx_Patch1(1)*dGdy_Patch2(2)-dGdy_Patch2(1)*dGdx_Patch1(2));
c=dGdx_Patch1(1)*dGdy_Patch1(2)-dGdy_Patch1(1)*dGdx_Patch1(2);




end

