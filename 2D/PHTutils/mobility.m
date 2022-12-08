function [ M,dM,d2M] = mobility( c )
%compute mobility for Cahn Hilliard equation

M=c.*(1-c);
dM=1-2.*c;
d2M=-2;
end

