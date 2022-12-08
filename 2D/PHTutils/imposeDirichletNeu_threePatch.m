function [ stiff, rhs, bcdof, bcval] = imposeDirichletNeu_threePatch(stiff, rhs, PHUTelem,GIFTmesh, p, q,P,type2Basis)
% impose neuman boundary condition for 5 patc hc1 example

%define side node indices
down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);

type2Basis_x=2*type2Basis-1;
type2Basis_y=2*type2Basis;

bcdof_left=[];
for i=1:length(PHUTelem{3})
    if isempty(PHUTelem{3}(i).children) && isempty(PHUTelem{3}(i).neighbor_left)
        
        bcdof_left = [bcdof_left,PHUTelem{3}(i).nodesGlobal(left_nodes)];
    end
end

%bcdof_left

%bcdof_left = [];
bcdof_left = unique(bcdof_left);
bcdof_left_xy = [2*bcdof_left-1, 2*bcdof_left];

bcdof = [bcdof_left_xy];
bcdof = setdiff(bcdof,[type2Basis_x,type2Basis_y]);

bcval = zeros(size(bcdof));

neumann_right = [1,2;2,2];
neumann = neumann_right;
num_neumann_elem = size(neumann,1);

%Gauss points
ngauss_edge = p+2;
[gauss_weight_edge, gauss_coord_edge] = quadrature( ngauss_edge, 'GAUSS', 1 );

gauss_coord_edge = gauss_coord_edge';

for i_neu=1:num_neumann_elem
    
    patchIndex = neumann(i_neu,1);
    
    elemList = sortEdgeElem( PHUTelem{patchIndex}, neumann(i_neu, 2));
    elemList = cell2mat(elemList);
    for i=1:length(elemList) %elem index
  
        corient = neumann(i_neu, 2);
        
        xmin = PHUTelem{patchIndex}(elemList(i)).vertex(1);
        xmax = PHUTelem{patchIndex}(elemList(i)).vertex(3);
        ymin = PHUTelem{patchIndex}(elemList(i)).vertex(2);
        ymax = PHUTelem{patchIndex}(elemList(i)).vertex(4);
        
        if (corient == 1) || (corient == 3)
            scalefac = (xmax-xmin)/2;
        else
            scalefac = (ymax-ymin)/2;
        end
        
        nument = size(PHUTelem{patchIndex}(elemList(i)).modifiedC,1);
        scrtx = PHUTelem{patchIndex}(elemList(i)).nodesGlobal(1:nument);
        
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
            
            [coord, dxdxi] = paramMap( GIFTmesh{patchIndex}, u_hat, v_hat, xmin, ymin, xmax, ymax);
            
            %[ coord,~,~, dxdxi] = paramMapPHUT( u_hat, v_hat, verts);
            plot(coord(1),coord(2),'+k')
            hold on
            drawnow
            
            
            %evaluate the basis functions
            [R, dR] = phtBasis(u_hat, v_hat, PHUTelem{patchIndex}(elemList(i)).modifiedC, p, q);
            
            
            %  jacobian of edge mapping;
            if((corient==1)||(corient==3))
                J = hypot(dxdxi(1,1), dxdxi(2,1));
            else
                J = hypot(dxdxi(1,2), dxdxi(2,2));
            end
            
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
            
            tmp = hypot(nor(1),nor(2));
            normal = nor/tmp; % normal vector in two dimensions
            
            taux = P;
            tauy = 0;
            localrhsed(1:2:end-1) = localrhsed(1:2:end-1) + R'.*taux.*scalefac.*gauss_weight_edge(igauss).*J;
            localrhsed(2:2:end) = localrhsed(2:2:end) + R'.*tauy.*scalefac.*gauss_weight_edge(igauss).*J;
            
        end
        rhs(dscrtx)=rhs(dscrtx)+localrhsed;
    end
end

%update the stiffness matrix and rhs with Dirichlet boundary values
[stiff,rhs] = feaplyc2sym(stiff,rhs, bcdof,bcval);

end

