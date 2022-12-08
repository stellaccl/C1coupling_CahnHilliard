function [dx,dy,d4x,d4y,d2xd2y]=computeDerivativesAnnalusBiharmonic(x,y)
%function []=computeDerivativesAnnalusBiharmonic()
% Rin=0.2;
% Rout=1;
% 
% syms x y
% % 
% u=(x^2+y^2-Rin^2)^2*(x^2+y^2-Rout^2)^2;

f_dx =@(x,y) 4*x*(x^2 + y^2 - 1)^2*(x^2 + y^2 - 1/25) + 4*x*(x^2 + y^2 - 1)*(x^2 + y^2 - 1/25)^2;
f_dy =@(x,y) 4*y*(x^2 + y^2 - 1)^2*(x^2 + y^2 - 1/25) + 4*y*(x^2 + y^2 - 1)*(x^2 + y^2 - 1/25)^2;

f_d4x=@(x,y) 576*x^2*(x^2 + y^2 - 1) + 576*x^2*(x^2 + y^2 - 1/25) + 24*(x^2 + y^2 - 1)^2 + 24*(x^2 + y^2 - 1/25)^2 + 96*(x^2 + y^2 - 1)*(x^2 + y^2 - 1/25) + 384*x^4;
f_d4y=@(x,y) 576*y^2*(x^2 + y^2 - 1) + 576*y^2*(x^2 + y^2 - 1/25) + 24*(x^2 + y^2 - 1)^2 + 24*(x^2 + y^2 - 1/25)^2 + 96*(x^2 + y^2 - 1)*(x^2 + y^2 - 1/25) + 384*y^4;

f_d2xd2y=@(x,y) 384*x^2*y^2 + 96*x^2*(x^2 + y^2 - 1) + 96*x^2*(x^2 + y^2 - 1/25) + 96*y^2*(x^2 + y^2 - 1) + 96*y^2*(x^2 + y^2 - 1/25) + 8*(x^2 + y^2 - 1)^2 + 8*(x^2 + y^2 - 1/25)^2 + 32*(x^2 + y^2 - 1)*(x^2 + y^2 - 1/25);

dx=f_dx(x,y);
dy=f_dy(x,y);
d4x=f_d4x(x,y);
d4y=f_d4y(x,y);
d2xd2y=f_d2xd2y(x,y);

end

