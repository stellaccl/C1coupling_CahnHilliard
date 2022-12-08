function [ stiff, rhs, bcdof, bcval ] = imposeDirichletNeuPHMPPlate2Patch_c1_degElev(stiff, rhs, PHTelem, p, q, numPatches,type2Basis)
%impose Dirichlet boundary conditions for the simply supported plate problem
%fixed (homogeneous) boundary conditions on the left, right, top, bottom sides
%for one patch only
% 2 nodes at each boundary
% two patches

%detect which nodes have support on the left and right boundary
bcdof_left = [];
bcdof_right = [];
bcdof_up = [];
bcdof_down = [];

%define side node indices
down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);

down_nodes2 = p+2:p+(p+2);
right_nodes2 = (p):(p+1):(p+1)*(q+1)-1;
up_nodes2 = (p+1)*q-p:(p+1)*(q+1)-p-1;
left_nodes2 = 2:(p+1):(1+(p+1)*q)+1;

%set the boundary degree of freedom and elements from the 1st patch
for indexPatch = 1:numPatches
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            
            numBasis=length(PHTelem{indexPatch}(i).nodesGlobal);
            
            if isempty(PHTelem{indexPatch}(i).neighbor_right) && (indexPatch==2)
                nodes=[right_nodes,right_nodes2];
                %                 nodes=[right_nodes];
                for ii=1:numBasis
                    if ~ismember(PHTelem{indexPatch}(i).nodesGlobal(ii),type2Basis)
                        nodesSupport=PHTelem{indexPatch}(i).modifiedC(ii,nodes);
                        if norm(nodesSupport)~=0
                            bcdof_right = [bcdof_right, PHTelem{indexPatch}(i).nodesGlobal(ii)];
                        end
                    end
                end
                %                 nodes=[right_nodes2];
                %                 for ii=1:numBasis
                %                     if ~ismember(PHTelem{indexPatch}(i).nodesGlobal(ii),type2Basis)
                %                         nodesSupport=PHTelem{indexPatch}(i).modifiedC(ii,nodes);
                %                         if norm(nodesSupport)~=0
                %                             bcdof_right = [bcdof_right, PHTelem{indexPatch}(i).nodesGlobal(ii)];
                %                         end
                %                     end
                %                 end
                
            end
            
            if isempty(PHTelem{indexPatch}(i).neighbor_down)
                nodes=[down_nodes,down_nodes2];
                %                 nodes=[down_nodes];
                for ii=1:numBasis
                    if ~ismember(PHTelem{indexPatch}(i).nodesGlobal(ii),type2Basis)
                        nodesSupport=PHTelem{indexPatch}(i).modifiedC(ii,nodes);
                        if norm(nodesSupport)~=0
                            bcdof_down = [bcdof_down, PHTelem{indexPatch}(i).nodesGlobal(ii)];
                        end
                    end
                end
                %                 nodes=[down_nodes2];
                %                 for ii=1:numBasis
                %                     if ~ismember(PHTelem{indexPatch}(i).nodesGlobal(ii),type2Basis)
                %                         nodesSupport=PHTelem{indexPatch}(i).modifiedC(ii,nodes);
                %                         if norm(nodesSupport)~=0
                %                             bcdof_down = [bcdof_down, PHTelem{indexPatch}(i).nodesGlobal(ii)];
                %                         end
                %                     end
                %                 end
            end
            
            if isempty(PHTelem{indexPatch}(i).neighbor_up)
                
                nodes=[up_nodes,up_nodes2];
                %                 nodes=[up_nodes];
                for ii=1:numBasis
                    if ~ismember(PHTelem{indexPatch}(i).nodesGlobal(ii),type2Basis)
                        nodesSupport=PHTelem{indexPatch}(i).modifiedC(ii,nodes);
                        if norm(nodesSupport)~=0
                            bcdof_up = [bcdof_up, PHTelem{indexPatch}(i).nodesGlobal(ii)];
                        end
                    end
                end
                
                %                 nodes=[up_nodes2];
                %                 for ii=1:numBasis
                %                     if ~ismember(PHTelem{indexPatch}(i).nodesGlobal(ii),type2Basis)
                %                         nodesSupport=PHTelem{indexPatch}(i).modifiedC(ii,nodes);
                %                         if norm(nodesSupport)~=0
                %                             bcdof_up = [bcdof_up, PHTelem{indexPatch}(i).nodesGlobal(ii)];
                %                         end
                %                     end
                %                 end
                
            end
            
            if isempty(PHTelem{indexPatch}(i).neighbor_left) && (indexPatch==1)
                
                nodes=[left_nodes,left_nodes2];
                %                 nodes=[left_nodes];
                for ii=1:numBasis
                    if ~ismember(PHTelem{indexPatch}(i).nodesGlobal(ii),type2Basis)
                        nodesSupport=PHTelem{indexPatch}(i).modifiedC(ii,nodes);
                        if norm(nodesSupport)~=0
                            bcdof_left = [bcdof_left, PHTelem{indexPatch}(i).nodesGlobal(ii)];
                        end
                    end
                end
                %                 nodes=[left_nodes2];
                %                 for ii=1:numBasis
                %                     if ~ismember(PHTelem{indexPatch}(i).nodesGlobal(ii),type2Basis)
                %                         nodesSupport=PHTelem{indexPatch}(i).modifiedC(ii,nodes);
                %                         if norm(nodesSupport)~=0
                %                             bcdof_left = [bcdof_left, PHTelem{indexPatch}(i).nodesGlobal(ii)];
                %                         end
                %                     end
                %                 end
                
            end
        end
    end
end


%remove duplicated entries
bcdof = unique([bcdof_down, bcdof_up, bcdof_left, bcdof_right]);
%make sure the type2basis are not included in the Dirichlet boundary
%conditions
bcdof = setdiff(bcdof,type2Basis);

bcval = zeros(size(bcdof));

%update the stiffness matrix and rhs with Dirichlet boundary values
%bcdof
%bcval
[stiff,rhs]=feaplyc2sym(stiff, rhs, bcdof, bcval);
