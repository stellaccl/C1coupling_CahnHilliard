function [ coord, dxdxi, d2xdxi2, dxdxi2] = paramMapShell_multiElemGIFTmesh( indexElem,GIFTmesh, u_hat, v_hat, xmin, ymin, xmax, ymax)
%maps (u_hat, v_hat) from the reference coordinate, to physical space using the geometry
%information stored in GIFTmesh

%map (u_hat, v_hat) on [-1, 1]x[-1,1] to (xi, eta) on [xmin, xmax]x[ymin,
%ymax]

xi = u_hat*(xmax-xmin)/2+(xmax+xmin)/2;
eta = v_hat*(ymax-ymin)/2+(ymax+ymin)/2;


%locate the vertices of the element on the coarse mesh
for inE = 1:GIFTmesh.numberElements

    if (xi<=GIFTmesh.elementVertex(inE,3)) && (xi>=GIFTmesh.elementVertex(inE,1)) && (eta<=GIFTmesh.elementVertex(inE,4)) && (eta>=GIFTmesh.elementVertex(inE,2))

        %disp('enter')

        e = inE;
        gxmin = GIFTmesh.elementVertex(e,1);
        gxmax = GIFTmesh.elementVertex(e,3);
        gymin = GIFTmesh.elementVertex(e,2);
        gymax = GIFTmesh.elementVertex(e,4);

    end
end

% e = indexElem;
% gxmin = GIFTmesh.elementVertex(e,1);
% gxmax = GIFTmesh.elementVertex(e,3);
% gymin = GIFTmesh.elementVertex(e,2);
% gymax = GIFTmesh.elementVertex(e,4);

%map (xi, eta) on [gxmin, gxmax]x[gymin, gymax] to (uu_hat, vv_hat) on [-1,
%1];

uu_hat = (2*xi - gxmin - gxmax)/(gxmax-gxmin);
vv_hat = (2*eta - gymin - gymax)/(gymax-gymin);

nodes = GIFTmesh.elementNode(e,:);

wgts = GIFTmesh.c_net(nodes, 4);
[N, dN, ddN] = nurbshape2d_gift_plate(uu_hat, vv_hat, wgts, GIFTmesh.C(:,:,e), GIFTmesh.p, GIFTmesh.q);
cpts = GIFTmesh.c_net(nodes,1:3);

%multiply by the jacobian of the transformation from reference
%space to the parameter space
dN(1,:) = dN(1,:)*2/(gxmax-gxmin);
dN(2,:) = dN(2,:)*2/(gymax-gymin);

ddN(1,:) = ddN(1,:)*(2/(gxmax-gxmin))^2;
ddN(2,:) = ddN(2,:)*(2/(gxmax-gxmin))*(2/(gymax-gymin));
ddN(3,:) = ddN(3,:)*(2/(gymax-gymin))^2;

%calculate the coordinates in the physical space
coord = N*cpts;

%calculate the Jacobian of the transformation


dxdxi = dN*cpts;


% Set up the second derivatives matrix and the matrix of squared first derivatives
d2xdxi2 = ddN*cpts;

dxdxi2 = [dxdxi(1,1)^2, 2*dxdxi(1,1)*dxdxi(1,2), dxdxi(1,2)^2;...
    dxdxi(1,1)*dxdxi(2,1), dxdxi(1,1)*dxdxi(2,2)+dxdxi(1,2)*dxdxi(2,1), dxdxi(1,2)*dxdxi(2,2);...
    dxdxi(2,1)^2, 2*dxdxi(2,1)*dxdxi(2,2), dxdxi(2,2)^2];





