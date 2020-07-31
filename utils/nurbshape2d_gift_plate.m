function [N, dN, ddN]=nurbshape2d_gift_plate(u_hat, v_hat, wgts, Ce, p, q)

%calculate the shape function and up to second derivatives
%INPUT: e - element index
%       u_hat - evaluation point in u-direction reference coordinates (from -1 to 1)
%       v_h - evaluation point in v-direction reference coordinates (from -1 to 1)
%       Ce - local Bezier extraction operator
%       p,q - polynomial degrees
%OUTPUT: N - value of the p+1 shape functions at the evaluation point
%        dN - derivatives of the p+1 shape functions in physical
%        coordinates (with respect to x)
%       ddN - 2nd derivatives of the p+1 shape functions in physical
%       coordinates (with respect to x)
% ------------------------------------------------------------------

%     evaluate 1d shape functions and derivatives 
[B_u, dB_u, ddB_u] = bernstein_basis(u_hat,p);
[B_v, dB_v, ddB_v] = bernstein_basis(v_hat,q);

B = zeros(1, (p+1)*(q+1));
dBxi = zeros(1, (p+1)*(q+1));
dBeta = zeros(1, (p+1)*(q+1));
ddBxi = zeros(1, (p+1)*(q+1));
ddBeta = zeros(1, (p+1)*(q+1));
ddBxieta = zeros(1, (p+1)*(q+1));

basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        
        %  Bernstein basis functions;    
        B(basisCounter) = B_u(i)*B_v(j);
        dBxi(basisCounter) = dB_u(i)*B_v(j);
        dBeta(basisCounter) = B_u(i)*dB_v(j);        
        ddBxi(basisCounter) = ddB_u(i)*B_v(j);
        ddBeta(basisCounter) = B_u(i)*ddB_v(j);        
        ddBxieta(basisCounter) = dB_u(i)*dB_v(j);        
    end
end

%form B-spline basis functions using the Bezier extraction operator
B = B*Ce';
dBxi = dBxi*Ce';
dBeta = dBeta*Ce';
ddBxi = ddBxi*Ce';
ddBeta = ddBeta*Ce';
ddBxieta = ddBxieta*Ce';

N = zeros(1, (p+1)*(q+1));
dN = zeros(2, (p+1)*(q+1));
ddN = zeros(3, (p+1)*(q+1));

% Multiply each B-spline function with corresponding weight
N(1,:) = B .* wgts';
dN(1,:) = dBxi .* wgts';
dN(2,:) = dBeta .* wgts';
ddN(1,:) = ddBxi .* wgts';
ddN(2,:) = ddBxieta .* wgts';
ddN(3,:) = ddBeta .* wgts';

% Compute the sums of B-spline functions
w_sum = sum(N(1,:));
dw_xi = sum(dN(1,:));
dw_eta = sum(dN(2,:));
d2w_xi = sum(ddN(1,:));
d2w_xieta = sum(ddN(2,:));
d2w_eta = sum(ddN(3,:));

% Compute NURBS basis functions and its first and second derivatives in
% local coordinates
ddN(1,:) = ddN(1,:)/w_sum - (2*dN(1,:)*dw_xi + N*d2w_xi)/w_sum^2 + 2*N*dw_xi^2/w_sum^3; %dxidxi derivative
ddN(3,:) = ddN(3,:)/w_sum - (2*dN(2,:)*dw_eta + N*d2w_eta)/w_sum^2 + 2*N*dw_eta^2/w_sum^3; %detadeta derivative
ddN(2,:) = ddN(2,:)/w_sum - (dN(1,:)*dw_eta + dN(2,:)*dw_xi + N*d2w_xieta)/w_sum^2 + ...
    2*N*dw_xi*dw_eta/w_sum^3; %dxideta derivative

dN(1,:) = dN(1,:)/w_sum - N*dw_xi/w_sum^2;
dN(2,:) = dN(2,:)/w_sum - N*dw_eta/w_sum^2;

N = N/w_sum;

