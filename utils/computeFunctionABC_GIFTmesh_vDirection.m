function [a,b,c] = computeFunctionABC_GIFTmesh_vDirection( PHUTelem,GIFTmesh,indexPatch1,indexElem1,indexPatch2,indexElem2,vref)

xmin = PHUTelem{indexPatch1}(indexElem1).vertex(1);
xmax = PHUTelem{indexPatch1}(indexElem1).vertex(3);
ymin = PHUTelem{indexPatch1}(indexElem1).vertex(2);
ymax = PHUTelem{indexPatch1}(indexElem1).vertex(4);
uref=1;
[coord, dxdxi] = paramMap_rotated( GIFTmesh{indexPatch1}, uref, vref, xmin, ymin, xmax, ymax);
% plot(coord(1),coord(2),'+r')
% hold on
%[~, dxdxi] = paramMap( GIFTmesh{indexPatch1}, uref, vref, xmin, ymin, xmax, ymax);
dGdx_Patch1Elem1=dxdxi(1,:);
dGdy_Patch1Elem1=dxdxi(2,:);

xmin = PHUTelem{indexPatch2}(indexElem2).vertex(1);
xmax = PHUTelem{indexPatch2}(indexElem2).vertex(3);
ymin = PHUTelem{indexPatch2}(indexElem2).vertex(2);
ymax = PHUTelem{indexPatch2}(indexElem2).vertex(4);
uref=-1;
[coord, dxdxi] = paramMap_rotated( GIFTmesh{indexPatch2}, uref, vref, xmin, ymin, xmax, ymax);
% plot(coord(1),coord(2),'+b')
% hold on
%[~, dxdxi] = paramMap( GIFTmesh{indexPatch2}, uref, vref, xmin, ymin, xmax, ymax);
dGdx_Patch2Elem1=dxdxi(1,:);
dGdy_Patch2Elem1=dxdxi(2,:);

a=dGdx_Patch1Elem1(1)*dGdy_Patch2Elem1(2)-dGdy_Patch2Elem1(1)*dGdx_Patch1Elem1(2);
b=dGdx_Patch1Elem1(1)*dGdx_Patch2Elem1(2)-dGdx_Patch2Elem1(1)*dGdx_Patch1Elem1(2);
c=dGdx_Patch2Elem1(1)*dGdy_Patch2Elem1(2)-dGdy_Patch2Elem1(1)*dGdx_Patch2Elem1(2);

% dGdx_Patch1Elem1
% dGdy_Patch1Elem1
% 
% dGdx_Patch2Elem1
% dGdy_Patch2Elem1
% pause

end

