function [ coord, dxdxi,d2xdxi2,dxdxi2] = paramMap1D( GIFTmesh, u_hat, xmin,xmax)

%maps (u_hat, v_hat) from the reference coordinate, to physical space using the geometry
%information stored in GIFTmesh

%map (u_hat, v_hat) on [-1, 1]x[-1,1] to (xi, eta) on [xmin, xmax]
xi = u_hat*(xmax-xmin)/2+(xmax+xmin)/2;

gxmin = GIFTmesh.elementVertex(1);
gxmax = GIFTmesh.elementVertex(2);

%map (xi, eta) on [gxmin, gxmax]x[gymin, gymax] to (uu_hat, vv_hat) on [-1,
%1];
uu_hat = (2*xi - gxmin - gxmax)/(gxmax-gxmin);

nodes=GIFTmesh.elementNode ;

wgts = GIFTmesh.c_net(nodes,2);

[N, dN,ddN] = nurbshape1d_gift(uu_hat, wgts, GIFTmesh.C, GIFTmesh.p);

cpts = GIFTmesh.c_net(nodes,1);
%multiply by the jacobian of the transformation from reference
%space to the parameter space
dN(1,:) = dN(1,:)*2/(gxmax-gxmin);
%ddN(1,:) = ddN(1,:)*(2/(gxmax-gxmin))^2;

%calculate the coordinates in the physical space
coord = N*cpts;

%calculate the Jacobian of the transformation
dxdxi = dN*cpts;

% Set up the second derivatives matrix and the matrix of squared first derivatives
d2xdxi2 = ddN*cpts;
% % 
dxdxi2 = dxdxi^2;



