function [ d4x ] = computeDerivatives1DBiharmonic(x)

 %syms x 
 %u=-(x*(1-x))^2;
d4x=-24;
 %diff(u,x,4)
end

