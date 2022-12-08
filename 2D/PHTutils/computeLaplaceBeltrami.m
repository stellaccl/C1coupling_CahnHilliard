%function [ laplacian] = computeLaplaceBeltrami(G,dxdu,dxdv,dydu,dydv,dzdu,dzdv,d2xdu2,d2xdv2,d2xdvdu,d2ydu2,d2ydv2,d2ydvdu,d2zdu2,d2zdv2,d2zdvdu,dphidu,dphidv,d2phidu2,d2phidv2,d2phidvdu)
function [ laplacian] = computeLaplaceBeltrami(G,dxdxi,d2xdxi2,dphidu,dphidv,d2phidu2,d2phidv2,d2phidvdu)
%compute laplacian according to "On Isogeometric Subdivision Methods for
%PDEs on Surfaces by Bert Juttler"

dxdu=dxdxi(1,1);
dxdv=dxdxi(2,1);
dydu=dxdxi(1,2);
dydv=dxdxi(2,2);
dzdu=dxdxi(1,3);
dzdv=dxdxi(2,3);

d2xdu2=d2xdxi2(1,1);
d2ydu2=d2xdxi2(1,2);
d2zdu2=d2xdxi2(1,3);

d2xdvdu=d2xdxi2(2,1);
d2ydvdu=d2xdxi2(2,2);
d2zdvdu=d2xdxi2(2,3);

d2xdv2=d2xdxi2(3,1);
d2ydv2=d2xdxi2(3,2);
d2zdv2=d2xdxi2(3,3);

part1=-(1/2)*(((dxdv)^2+(dydv)^2+(dzdv)^2)*(dphidu)+(-(dxdu)*(dxdv)-(dydu)*(dydv)-(dzdu)*(dzdv))*(dphidv))*((2*(dxdu)*(d2xdu2)+2*(dydu)*(d2ydu2)+2*(dzdu)*(d2zdu2))*((dxdv)^2+(dydv)^2+(dzdv)^2)+((dxdu)^2+(dydu)^2+(dzdu)^2)*(2*(dxdv)*(d2xdvdu)+2*(dydv)*(d2ydvdu)+2*(dzdv)*(d2zdvdu))-(2*((dxdu)*(dxdv)+(dydu)*(dydv)+(dzdu)*(dzdv)))*((d2xdu2)*(dxdv)+(dxdu)*(d2xdvdu)+(d2ydu2)*(dydv)+(dydu)*(d2ydvdu)+(d2zdu2)*(dzdv)+(dzdu)*(d2zdvdu)))/(((dxdu)^2+(dydu)^2+(dzdu)^2)*((dxdv)^2+(dydv)^2+(dzdv)^2)-((dxdu)*(dxdv)+(dydu)*(dydv)+(dzdu)*(dzdv))^2)^(3/2)+((2*(dxdv)*(d2xdvdu)+2*(dydv)*(d2ydvdu)+2*(dzdv)*(d2zdvdu))*(dphidu)+((dxdv)^2+(dydv)^2+(dzdv)^2)*(d2phidu2)+(-(d2xdu2)*(dxdv)-(dxdu)*(d2xdvdu)-(d2ydu2)*(dydv)-(dydu)*(d2ydvdu)-(d2zdu2)*(dzdv)-(dzdu)*(d2zdvdu))*(dphidv)+(-(dxdu)*(dxdv)-(dydu)*(dydv)-(dzdu)*(dzdv))*(d2phidvdu))/sqrt(((dxdu)^2+(dydu)^2+(dzdu)^2)*((dxdv)^2+(dydv)^2+(dzdv)^2)-((dxdu)*(dxdv)+(dydu)*(dydv)+(dzdu)*(dzdv))^2);
part2=-(1/2)*((-(dxdu)*(dxdv)-(dydu)*(dydv)-(dzdu)*(dzdv))*(dphidu)+((dxdu)^2+(dydu)^2+(dzdu)^2)*(dphidv))*((2*(dxdu)*(d2xdvdu)+2*(dydu)*(d2ydvdu)+2*(dzdu)*(d2zdvdu))*((dxdv)^2+(dydv)^2+(dzdv)^2)+((dxdu)^2+(dydu)^2+(dzdu)^2)*(2*(dxdv)*(d2xdv2)+2*(dydv)*(d2ydv2)+2*(dzdv)*(d2zdv2))-(2*((dxdu)*(dxdv)+(dydu)*(dydv)+(dzdu)*(dzdv)))*((dxdv)*(d2xdvdu)+(dxdu)*(d2xdv2)+(dydv)*(d2ydvdu)+(dydu)*(d2ydv2)+(dzdv)*(d2zdvdu)+(dzdu)*(d2zdv2)))/(((dxdu)^2+(dydu)^2+(dzdu)^2)*((dxdv)^2+(dydv)^2+(dzdv)^2)-((dxdu)*(dxdv)+(dydu)*(dydv)+(dzdu)*(dzdv))^2)^(3/2)+((-(dxdv)*(d2xdvdu)-(dxdu)*(d2xdv2)-(dydv)*(d2ydvdu)-(dydu)*(d2ydv2)-(dzdv)*(d2zdvdu)-(dzdu)*(d2zdv2))*(dphidu)+(-(dxdu)*(dxdv)-(dydu)*(dydv)-(dzdu)*(dzdv))*(d2phidvdu)+(2*(dxdu)*(d2xdvdu)+2*(dydu)*(d2ydvdu)+2*(dzdu)*(d2zdvdu))*(dphidv)+((dxdu)^2+(dydu)^2+(dzdu)^2)*(d2phidv2))/sqrt(((dxdu)^2+(dydu)^2+(dzdu)^2)*((dxdv)^2+(dydv)^2+(dzdv)^2)-((dxdu)*(dxdv)+(dydu)*(dydv)+(dzdu)*(dzdv))^2);

laplacian=1/sqrt(det(G))*(part1+part2);

end
