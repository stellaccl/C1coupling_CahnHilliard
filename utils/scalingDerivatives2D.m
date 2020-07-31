function [dRdu_patch1,  dRdv_patch1,dRdu_patch2,  dRdv_patch2] =scalingDerivatives2D(PHUTelem,indexPatch1,indexPatch2,indexElem1,indexElem2,dRdu_patch1,  dRdv_patch1, dRdu_patch2,  dRdv_patch2)

%scale the derivatives with respect to the transformation
%from the reference element to the parameter space
xmin = PHUTelem{indexPatch1}(indexElem1).vertex(1);
xmax = PHUTelem{indexPatch1}(indexElem1).vertex(3);
ymin = PHUTelem{indexPatch1}(indexElem1).vertex(2);
ymax = PHUTelem{indexPatch1}(indexElem1).vertex(4);

dRdu_patch1 = dRdu_patch1*2/(xmax-xmin);
dRdv_patch1 = dRdv_patch1*2/(ymax-ymin);

xmin = PHUTelem{indexPatch2}(indexElem2).vertex(1);
xmax = PHUTelem{indexPatch2}(indexElem2).vertex(3);
ymin = PHUTelem{indexPatch2}(indexElem2).vertex(2);
ymax = PHUTelem{indexPatch2}(indexElem2).vertex(4);

dRdu_patch2= dRdu_patch2*2/(xmax-xmin);
dRdv_patch2= dRdv_patch2*2/(ymax-ymin);

end

