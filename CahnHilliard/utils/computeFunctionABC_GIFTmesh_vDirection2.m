function [a,b,c] = computeFunctionABC_GIFTmesh_vDirection2( PHUTelemA,PHUTelemB,GIFTmeshA,GIFTmeshB,vref,urefA,uferB)
%use for 5patches example 

xmin = PHUTelemA.vertex(1);
xmax = PHUTelemA.vertex(3);
ymin = PHUTelemA.vertex(2);
ymax = PHUTelemA.vertex(4);
uref=1;

[coord,dxdxi] = paramMap(GIFTmeshA,uref,vref,xmin,ymin,xmax,ymax);
             
dGdx_Patch1Elem1=dxdxi(1,:);
dGdy_Patch1Elem1=dxdxi(2,:);

xmin = PHUTelemB.vertex(1);
xmax = PHUTelemB.vertex(3);
ymin = PHUTelemB.vertex(2);
ymax = PHUTelemB.vertex(4);
uref=-1;
[coord,dxdxi] = paramMap( GIFTmeshB,uref,vref,xmin,ymin,xmax,ymax);

dGdx_Patch2Elem1=dxdxi(1,:);
dGdy_Patch2Elem1=dxdxi(2,:);

a=dGdx_Patch1Elem1(1)*dGdy_Patch2Elem1(2)-dGdy_Patch2Elem1(1)*dGdx_Patch1Elem1(2);
b=dGdx_Patch1Elem1(1)*dGdx_Patch2Elem1(2)-dGdx_Patch2Elem1(1)*dGdx_Patch1Elem1(2);
c=dGdx_Patch2Elem1(1)*dGdy_Patch2Elem1(2)-dGdy_Patch2Elem1(1)*dGdx_Patch2Elem1(2);

end

