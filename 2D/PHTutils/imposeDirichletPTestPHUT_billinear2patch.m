function [ stiff, rhs, bcdof, bcval] = imposeDirichletPTestPHUT_billinear2patch(stiff, rhs, PHUTelem, p, q, bound_disp,patchBoundaries)
%impose Dirichlet boundary conditions for the patch test problem for PHUT
%splines
%(elasticity equation)
%fixed (homogeneous) on the left side, bound_disp displacement on the right
%for bilinear example

%define side node indices
down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);

% bcdof_left = PHUTelem(3).nodes(left_nodes);
% bcdof_right = [PHUTelem(1).nodes(right_nodes), PHUTelem(2).nodes(right_nodes)];
numBoundaries = size(patchBoundaries,1);
bcdof_left=[];
bcdof_right=[];

for boundaryIndex = 1:numBoundaries
    patchA = patchBoundaries(boundaryIndex,1);
    patchB = patchBoundaries(boundaryIndex,2);
    edgeA = 4;
    edgeB = 2;

    patchA=cell2mat(patchA);
    patchB=cell2mat(patchB);
    
    [elementsA] = sortEdgeElem( PHUTelem{patchA}, edgeA);
    [elementsB] = sortEdgeElem( PHUTelem{patchB}, edgeB);
    
    elementsA=elementsA{1};
    elementsB=elementsB{1};

    for i=1:length(elementsA)
        bcdof_left = [bcdof_left, PHUTelem{1}(elementsA(i)).nodesGlobal(left_nodes)];
    end
    
    for i=1:length(elementsB)
        bcdof_right = [bcdof_right, PHUTelem{2}(elementsB(i)).nodesGlobal(right_nodes)];
    end
end


% bcdof_left = [ PHUTelem{1}(2).nodesGlobal(left_nodes),PHUTelem{1}(4).nodesGlobal(left_nodes)];
% bcdof_right = [ PHUTelem{2}(3).nodesGlobal(right_nodes),PHUTelem{2}(5).nodesGlobal(right_nodes)];
bcdof_right = unique(bcdof_right);
bcdof_left = unique(bcdof_left);

%impose the constraints (zero x and y displacement) in the linear system
bcdof_left_xy = [2*bcdof_left-1, 2*bcdof_left];
bcdof_right_x = 2*bcdof_right-1;
bcdof_right_y = 2*bcdof_right;

bcval_left_xy = zeros(size(bcdof_left_xy));
bcval_right_x = bound_disp*ones(size(bcdof_right_x));
bcval_right_y = zeros(size(bcdof_right_y));

bcdof = [bcdof_left_xy, bcdof_right_x, bcdof_right_y];
bcval = [bcval_left_xy, bcval_right_x, bcval_right_y];

%update the stiffness matrix and rhs with Dirichlet boundary values
[stiff,rhs] = feaplyc2sym(stiff,rhs, bcdof,bcval);

