function  [ stiff, rhs, bcdof, bcval] = imposeDirichletBiharmonic_1D(stiff, rhs, PHTelem, p,type2Basis)

%detect which nodes have support on the left and right boundary
bcdof_left = [];
bcdof_right = [];

%define side node indices
left_nodes = 1:2;
right_nodes = p:p+1;

%set the boundary degree of freedom
for indexPatch = 1:length(PHTelem)
    for i=1:length(PHTelem{indexPatch})
        
        if indexPatch == 1 && i==1
             bcdof_left = [bcdof_left, PHTelem{indexPatch}(i).nodesGlobal(left_nodes)];
        end

        if indexPatch == 2 && i==length(PHTelem{indexPatch})
             bcdof_right = [bcdof_right, PHTelem{indexPatch}(i).nodesGlobal(right_nodes)];
        end

    end
end

bcdof = unique([bcdof_left, bcdof_right]);
bcdof = setdiff(bcdof,type2Basis);
bcval = zeros(size(bcdof));
[stiff,rhs]=feaplyc2sym(stiff, rhs, bcdof, bcval);


end

