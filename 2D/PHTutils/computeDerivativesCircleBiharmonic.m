function [dx,dy,d4x,d4y,d2xd2y]=computeDerivativesCircleBiharmonic(x,y)

% r=0.5;
% syms x y
% u=(x^2+y^2-r^2)^2;
% diff(u,x,x,y,y)

f_dx =@(x,y) 4*x*(x^2 + y^2 - 1/4);
f_dy =@(x,y) 4*y*(x^2 + y^2 - 1/4);

f_d4x=@(x,y) 24;
f_d4y=@(x,y) 24;

f_d2xd2y=@(x,y) 8;

dx=f_dx(x,y);
dy=f_dy(x,y);
d4x=f_d4x(x,y);
d4y=f_d4y(x,y);
d2xd2y=f_d2xd2y(x,y);

end

