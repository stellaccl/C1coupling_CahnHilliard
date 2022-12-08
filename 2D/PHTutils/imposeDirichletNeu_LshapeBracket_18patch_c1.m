function [ stiff, rhs, bcdof, bcval ] = imposeDirichletNeu_LshapeBracket_18patch_c1(stiff, rhs, PHTelem, GIFTmesh, p, q, tauy,type2Basis)
%impose Dirichlet boundary conditions for Lshape bracket
%fix the last two holes
%Neumann (traction) boundary conditions on the last hole


ngauss_edge = q+2;
[gauss_weight_edge, gauss_coord_edge] = quadrature( ngauss_edge, 'GAUSS', 1 );

%detect which nodes have support on the left and right boundary
bcdof_left = [];
bcdof_right = [];
bcdof_up = [];
bcdof_down = [];

%for each neumann edge store the element index and orientation
%orientation: 1-down, 2-right, 3-up, 4-left
neumann_left = [];
neumann_right = [];
neumann_up = [];
neumann_down =[];

%define side node indices
down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);

bcdof = [];

for indexPatch=11:18
    for indexElem=1:length(PHTelem{indexPatch})
        if  isempty(PHTelem{indexPatch}(indexElem).neighbor_left)
            bcdof=[bcdof, PHTelem{indexPatch}(indexElem).nodesGlobal(left_nodes)];
        end
    end
end

for indexPatch=1:4
    for indexElem=1:length(PHTelem{indexPatch})
        if  isempty(PHTelem{indexPatch}(indexElem).neighbor_left)
            neumann_left= [neumann_left; indexElem, 4, indexPatch];
        end
    end
end

type2Basis_x=2*type2Basis-1;
type2Basis_y=2*type2Basis;

bcdof = [2*bcdof, 2*bcdof-1];
bcdof=unique(bcdof);
bcdof = setdiff(bcdof,[type2Basis_x,type2Basis_y]);
bcval = zeros(size(bcdof));

%impose Neumann boundary conditons
neumann = neumann_left;

for i_neu=1:size(neumann,1)
    
    i = neumann(i_neu,1);
    patchIndex = neumann(i_neu,3);
    corient = neumann(i_neu, 2);
    
    xmin = PHTelem{patchIndex}(i).vertex(1);
    xmax = PHTelem{patchIndex}(i).vertex(3);
    ymin = PHTelem{patchIndex}(i).vertex(2);
    ymax = PHTelem{patchIndex}(i).vertex(4);
    
    if (corient == 1) || (corient == 3)
        scalefac = (xmax-xmin)/2;
    else
        scalefac = (ymax-ymin)/2;
    end
    
    nument = size(PHTelem{patchIndex}(i).modifiedC,1);
    
    
    scrtx = PHTelem{patchIndex}(i).nodesGlobal(1:nument);    
    dscrtx = reshape([2*scrtx-1; 2*scrtx],1,2*nument);
    localrhsed = zeros(2*nument, 1);
    
    %loop over Gauss points and compute the integral
    for igauss = 1:ngauss_edge
        
        %find the evaluation point, taking into account the edge type
        if corient==1
            v_hat = -1;
            u_hat = gauss_coord_edge(igauss);
        elseif corient==2
            u_hat = 1;
            v_hat = gauss_coord_edge(igauss);
        elseif corient==3
            v_hat = 1;
            u_hat = gauss_coord_edge(igauss);
        else
            u_hat = -1;
            v_hat = gauss_coord_edge(igauss);
        end
        
        %evaluate the basis functions
        R = phtBasis(u_hat, v_hat, PHTelem{patchIndex}(i).modifiedC, p, q);
        
        %evaluate the derivatives of the mapping from parameter
        %space to physical space
        [coord, dxdxi] = paramMap( GIFTmesh{patchIndex}, u_hat, v_hat, xmin, ymin, xmax, ymax);
        
        plot(coord(1),coord(2),'+r')
        hold on
        drawnow
        
        
        dxdxi = dxdxi';
        %  jacobian of edge mapping;
        if((corient==1)||(corient==3))
            ee = dxdxi(1,1)^2+dxdxi(2,1)^2;
        else
            ee = dxdxi(1,2)^2+dxdxi(2,2)^2;
        end
        
        % Jacobian of face mapping
        J = sqrt(ee);%Length of edge
        %------------------------------------------------------------------------;
        %  computation of normal;
        if(corient==1)
            nor(1) = dxdxi(2,1);% dy/dxi
            nor(2) = -dxdxi(1,1);%dx/dxi
        elseif(corient==2)
            nor(1) = dxdxi(2,2);
            nor(2) = -dxdxi(1,2);
        elseif(corient==3)
            nor(1) = -dxdxi(2,1);
            nor(2) = dxdxi(1,1);
        else
            nor(1) = -dxdxi(2,2);
            nor(2) = dxdxi(1,2);
        end
        
        tmp = sqrt(nor(1)^2 + nor(2)^2);
        normal = nor/tmp; % normal vector in two dimensions
                
        taux = 0;
        
        localrhsed(1:2:end-1) = localrhsed(1:2:end-1) + R'.*taux.*scalefac.*gauss_weight_edge(igauss).*J;
        localrhsed(2:2:end) = localrhsed(2:2:end) + R'.*tauy.*scalefac.*gauss_weight_edge(igauss).*J;
        
    end
    rhs(dscrtx)=rhs(dscrtx)+localrhsed;
end

[stiff,rhs]=feaplyc2(stiff, rhs, bcdof, bcval);

end


