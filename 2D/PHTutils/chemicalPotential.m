function [ mu,dmu,d2mu ] = chemicalPotential( c, theta,alpha )
%chemicalPotential for Chan hiliard equation

mu=0.5./theta.*log(c./(1-c))+1-2.*c;
mu=mu*3*alpha;

dmu=0.5./theta.*1./(c.*(1-c))-2;
dmu=dmu*3*alpha;
%
d2mu=0.5./theta.*(2.*c-1)/(c.*c.*(1-c).*(1-c));
d2mu=d2mu*3*alpha;
end

