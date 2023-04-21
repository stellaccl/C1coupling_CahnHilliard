function [N, dN,ddN]=nurbshape1d_gift(u_hat, wgts, Ce, p)

%calculate the shape function and first derivatives
%INPUT: e - element index
%       u_hat - evaluation point in u-direction reference coordinates (from -1 to 1)
%       v_h - evaluation point in v-direction reference coordinates (from -1 to 1)
%       Ce - local Bezier extraction operator
%       p,q - polynomial degrees
%OUTPUT: N - value of the p+1 shape functions at the evaluation point
%        dN - derivatives of the p+1 shape functions in physical
%        coordinates (with respect to x)
% ------------------------------------------------------------------

%     evaluate 1d shape functions and derivatives
[B_u, dB_u, ddB_u] = bernstein_basis(u_hat,p);

B = zeros(1, (p+1));
dBxi = zeros(1, (p+1));
ddBxi = zeros(1, (p+1));

basisCounter = 0;

for i=1:p+1
    basisCounter = basisCounter + 1;
    
    %  Bernstein basis functions;
    B(basisCounter) = B_u(i);
    dBxi(basisCounter) = dB_u(i);
    ddBxi(basisCounter) = ddB_u(i);
end

%form B-spline basis functions using the Bezier extraction operator
B = B*Ce';
dBxi = dBxi*Ce';
ddBxi = ddBxi*Ce';

N = zeros(1, (p+1));
dN = zeros(1, (p+1));
ddN = zeros(1, (p+1));

% Multiply each B-spline function with corresponding weight
N(1,:) = B .* wgts';
dN(1,:) = dBxi .* wgts';
ddN(1,:) = ddBxi .* wgts';

% Compute the sums of B-spline functions
w_sum = sum(N(1,:));
dw_xi = sum(dN(1,:));
d2w_xi = sum(ddN(1,:));

% Compute NURBS basis functions and its first derivatives in
% local coordinates
ddN(1,:) = ddN(1,:)/w_sum - (2*dN(1,:)*dw_xi + N*d2w_xi)/w_sum^2 + 2*N*dw_xi^2/w_sum^3; %dxidxi derivative

dN(1,:) = dN(1,:)/w_sum - N*dw_xi/w_sum^2;

N = N/w_sum;
