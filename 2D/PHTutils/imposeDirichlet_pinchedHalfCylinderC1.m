function [ stiff, rhs, bcdof, bcval, forcedNode ] = imposeDirichlet_pinchedHalfCylinderC1(stiff, rhs, PHTelem, p, q,sizeBasis, type2Basis)
%impose Dirichlet boundary conditions for the pinchedCylinder problem

%   z
%   |
% (D)------- (C)
%   |      |
%   |    L |
%   |  R   |
% (A)------- (B) --->x
% rigid diaphram: AB: u_x, u_y = 0
%Symmetry conditions on BC,CD and AD

bcdof_AB = [];
bcdof_CB = [];
bcdof_CD = [];
bcdof_AD = [];

%define side node indices
down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);

numPatches=length(PHTelem);
forcedNode=[];
%set the boundary degree of freedom and elements from the 1st patch
for indexPatch = 1:numPatches
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            
            if isempty(PHTelem{indexPatch}(i).neighbor_down)
                if length(PHTelem{indexPatch}(i).nodesGlobal)<(p+1)*(q+1)
                    bcdof_AB=[bcdof_AB,PHTelem{indexPatch}(i).nodesGlobal([1,2])];
                else
                    bcdof_AB=[bcdof_AB,PHTelem{indexPatch}(i).nodesGlobal(down_nodes)];
                end
            end
            
            if isempty(PHTelem{indexPatch}(i).neighbor_up)
                if length(PHTelem{indexPatch}(i).nodesGlobal)<(p+1)*(q+1)
                    bcdof_CD=[bcdof_CD,PHTelem{indexPatch}(i).nodesGlobal([10,11])];
                else
                    bcdof_CD=[bcdof_CD,PHTelem{indexPatch}(i).nodesGlobal(up_nodes)];
                    
                end
            end
            
            if isempty(PHTelem{indexPatch}(i).neighbor_left) && indexPatch==1
                if length(PHTelem{indexPatch}(i).nodesGlobal)<(p+1)*(q+1)
                    disp('error')
                    pause
                else
                    bcdof_CB=[bcdof_CB,PHTelem{indexPatch}(i).nodesGlobal(left_nodes)];
                end
            end
            
            if isempty(PHTelem{indexPatch}(i).neighbor_right) && indexPatch==2
                if length(PHTelem{indexPatch}(i).nodesGlobal)<(p+1)*(q+1)
                    disp('error')
                    pause
                else
                    bcdof_AD=[bcdof_AD,PHTelem{indexPatch}(i).nodesGlobal(right_nodes)];
                end
            end
            
            xmin = PHTelem{indexPatch}(i).vertex(1);
            xmax = PHTelem{indexPatch}(i).vertex(3);
            ymin = PHTelem{indexPatch}(i).vertex(2);
            ymax = PHTelem{indexPatch}(i).vertex(4);
            
            NEcorner=[xmax,ymax];
            SEcorner=[xmax,ymin];
            SWcorner=[xmin,ymin];
            NWcorner=[xmin,ymax];
            
            if indexPatch==1
                centerPoint=[1,0.5];
            elseif indexPatch==2
                centerPoint=[0,0.5];
            end
            
            if NEcorner==centerPoint
                forcedNode = [forcedNode,PHTelem{indexPatch}(i).nodesGlobal];
            elseif SEcorner==centerPoint
                forcedNode = [forcedNode,PHTelem{indexPatch}(i).nodesGlobal];
            elseif SWcorner==centerPoint
                forcedNode = [forcedNode,PHTelem{indexPatch}(i).nodesGlobal];
            elseif NWcorner==centerPoint
                forcedNode = [forcedNode,PHTelem{indexPatch}(i).nodesGlobal];
            end
            
            
        end
    end
end
forcedNode=unique(forcedNode);
xConsNodes = unique([bcdof_CD, bcdof_CB,bcdof_AD,bcdof_AB]);
yConsNodes = unique([bcdof_CD,bcdof_AB ]);
zConsNodes = [];

% xConsNodes = unique([bcdof_AB bcdof_AD]);
% yConsNodes = unique([bcdof_AB bcdof_CB]);
% zConsNodes = unique([bcdof_CD]);

udofs      = 3*xConsNodes-2; % global indecies  of the fixed x disps
vdofs      = 3*yConsNodes-1; % global indecies  of the fixed y disps
wdofs      = 3*zConsNodes;   % global indecies  of the fixed z disps

type2Basis_x=3*type2Basis-2;
type2Basis_y=3*type2Basis-1;
type2Basis_z=3*type2Basis;

%w     = 1e7;
%w=0;
%penaltyStiffness = w*[1 -1;-1 1];

for indexNodes=1:length(bcdof_CD)
    %sctr  = [bcdof_CD(indexNodes) bcdof_CD2(indexNodes)];
    sctr  = [bcdof_CD(indexNodes)];
    sctrx = 3*sctr-2;
    sctry = 3*sctr-1;
    sctrz = 3*sctr-0;
    
    stiff(sctrx,sctrx) = stiff(sctrx,sctrx);
    stiff(sctry,sctry) = stiff(sctry,sctry);
    stiff(sctrz,sctrz) = stiff(sctrz,sctrz);
    
end

for indexBasis=1:length(bcdof_CB)
    %sctr  = [bcdof_CB(indexBasis) bcdof_CB2(indexBasis)];
    sctr  = [bcdof_CB(indexBasis)];
    sctrx = 3*sctr-2;
    sctry = 3*sctr-1;
    sctrz = 3*sctr-0;
    
    stiff(sctrx,sctrx) = stiff(sctrx,sctrx);
    stiff(sctry,sctry) = stiff(sctry,sctry);
    stiff(sctrz,sctrz) = stiff(sctrz,sctrz);
end

for indexBasis=1:length(bcdof_AD)
    %sctr  = [bcdof_AD(indexBasis) bcdof_AD2(indexBasis)];
    sctr  = [bcdof_AD(indexBasis)];
    sctrx = 3*sctr-2;
    sctry = 3*sctr-1;
    sctrz = 3*sctr-0;
    
    stiff(sctrx,sctrx) = stiff(sctrx,sctrx);
    stiff(sctry,sctry) = stiff(sctry,sctry);
    stiff(sctrz,sctrz) = stiff(sctrz,sctrz);
end

%forcedNode = sizeBasis;
rhs(3*forcedNode-1)=-1/4;

%update the stiffness matrix and rhs with Dirichlet boundary values
%bcdof = [udofs, vdofs, wdofs];
bcdof = unique([udofs, vdofs, wdofs]);
bcdof = setdiff(bcdof,[type2Basis_x,type2Basis_y,type2Basis_z]);
bcval = zeros(size(bcdof));
[stiff,rhs]=feaplyc2sym(stiff, rhs, bcdof, bcval );
