function [ stiff, rhs, bcdof, bcval ] = imposeDirichletNeuPHMP_degElev(stiff, rhs, PHTelem, GIFTmesh, p, q, rad, tx)
%impose Dirichlet boundary conditions for elastic rectangle
%fixed (homogeneous) boundary conditions on the left side
%Neumann (traction) boundary conditions on the right side
%for use with the degree elevation subroutines

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

%set the boundary degree of freedom and elements from the 1st patch
for i=1:length(PHTelem{1})
    if isempty(PHTelem{1}(i).children)                
        if isempty(PHTelem{1}(i).neighbor_up)
            neumann_up = [neumann_up; i, 3, 1];
        end
        if isempty(PHTelem{1}(i).neighbor_left)
            bcdof_left = [bcdof_left, PHTelem{1}(i).nodesGlobal(left_nodes)];            
        end                
    end
end


%set the boundary degree of freedom and elements from the 2nd patch
for i=1:length(PHTelem{2})
    if isempty(PHTelem{2}(i).children)        
        if isempty(PHTelem{2}(i).neighbor_right)
            bcdof_right = [bcdof_right, PHTelem{2}(i).nodesGlobal(right_nodes)];  
        end
        if isempty(PHTelem{2}(i).neighbor_up)            
            neumann_up = [neumann_up; i, 3, 2];
        end        
    end
end


%remove duplicated entries
bcdof_right = unique(bcdof_right);
bcdof_left = unique(bcdof_left);

bcdof_left_y = 2*bcdof_left;

bcdof_right_x = 2*bcdof_right-1;

%impose the symmetry boundary conditions on the top and bottom edge
bcval_right_x = zeros(size(bcdof_right_x));
bcval_left_y = zeros(size(bcdof_left_y));

bcdof = [bcdof_right_x, bcdof_left_y];
bcval = [bcval_right_x, bcval_left_y];

%impose Neumann boundary conditons
neumann = neumann_up;
side_length = 0;
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
    
%     
%     if corient==1
%         scrtx = scrtx(down_nodes);
%         nument_edge = p+1;
%     elseif corient==2
%         scrtx = scrtx(right_nodes);
%         nument_edge = q+1;
%     elseif corient==3
%         scrtx = scrtx(up_nodes);
%         nument_edge = p+1;
%     else
%         scrtx = scrtx(left_nodes);
%         nument_edge = q+1;
%     end
%    dscrtx = reshape([2*scrtx-1; 2*scrtx],1,2*nument_edge);
    dscrtx = reshape([2*scrtx-1; 2*scrtx],1,2*nument);
    %localrhsed = zeros(2*nument_edge, 1);
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
        R = zeros(1,nument);
        for t=1:nument
            if PHTelem{patchIndex}(i).polyDegree(t)==p
                %evaluate the basis functions
                R(t) = phtBasis(u_hat, v_hat, PHTelem{patchIndex}(i).modifiedC(t,1:(p+1)*(q+1)), p, q);
            else
                R(t) = phtBasis(u_hat, v_hat, PHTelem{patchIndex}(i).modifiedC(t,:), p+1, q+1);
            end
        end
        
        %evaluate the derivatives of the mapping from parameter
        %space to physical space
        [coord, dxdxi] = paramMap( GIFTmesh{patchIndex}, u_hat, v_hat, xmin, ymin, xmax, ymax);
%         plot(coord(1),coord(2),'.k','MarkerSize',15)
%         hold on
       
        stress = ghole(coord(1), coord(2), rad, tx);
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
%         quiver(coord(1),coord(2),normal(1),normal(2))
%         hold on
%          drawnow

%         %consider on the values corresponding to the shape functions on the
%         %boundary
%         if corient==1
%             R = R(down_nodes)';
%         elseif corient==2
%             R = R(right_nodes)';
%         elseif corient==3
%             R = R(up_nodes)';
%         else
%             R = R(left_nodes)';
%         end
        
        taux = normal(1)*stress(1) + normal(2)*stress(3);
        tauy = normal(1)*stress(3) + normal(2)*stress(2);
 
        R = R';
        %R.*taux.*scalefac.*gauss_weight_edge(igauss).*J
        localrhsed(1:2:end-1) = localrhsed(1:2:end-1) + R.*taux.*scalefac.*gauss_weight_edge(igauss).*J;
        localrhsed(2:2:end) = localrhsed(2:2:end) + R.*tauy.*scalefac.*gauss_weight_edge(igauss).*J;
        side_length = side_length + scalefac.*gauss_weight_edge(igauss).*J;
       % pause
    end
    rhs(dscrtx)=rhs(dscrtx)+localrhsed;
end

%update the stiffness matrix and rhs with Dirichlet boundary values
[stiff,rhs]=feaplyc2(stiff, rhs, bcdof, bcval);
side_length
