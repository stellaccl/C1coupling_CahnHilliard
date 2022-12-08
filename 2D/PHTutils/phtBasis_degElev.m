function [N, dN]=phtBasis_degElev(PHTelem,patchIndex,i,u_hat, v_hat, p, q,nument)

%calculate the shape function and first derivatives
%INPUT: e - element index
%       u_hat - evaluation point in u-direction reference coordinates (from -1 to 1)
%       v_hat - evaluation point in v-direction reference coordinates (from -1 to 1)
%       Ce - local Bezier extraction operator
%       p,q - polynomial degrees
%OUTPUT: N - value of the p+1 shape functions at the evaluation point
%        dN - derivatives of the shape functions in %
% ------------------------------------------------------------------

%     evaluate 1d shape functions and derivatives
[B_u, dB_u] = bernstein_basis(u_hat,p);
[B_v, dB_v] = bernstein_basis(v_hat,q);

%     evaluate 1d shape functions and derivatives
[B_uDE, dB_uDE] = bernstein_basis(u_hat,p+1);
[B_vDE, dB_vDE] = bernstein_basis(v_hat,q+1);

B = zeros(1, (p+1)*(q+1));
dBxi = zeros(1, (p+1)*(q+1));
dBeta = zeros(1, (p+1)*(q+1));

B_DE = zeros(1, (p+2)*(q+2));
dBxi_DE = zeros(1, (p+2)*(q+2));
dBeta_DE = zeros(1, (p+2)*(q+2));

basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        
        %  Bernstein basis functions;
        B(basisCounter) = B_u(i)*B_v(j);
        dBxi(basisCounter) = dB_u(i)*B_v(j);
        dBeta(basisCounter) = B_u(i)*dB_v(j);
    end
end

basisCounter = 0;
for j=1:q+2
    for i=1:p+2
        basisCounter = basisCounter + 1;
        %  Bernstein basis functions;
        B_DE(basisCounter) = B_uDE(i)*B_vDE(j);
        dBxi_DE(basisCounter) = dB_uDE(i)*B_vDE(j);
        dBeta_DE(basisCounter) = B_uDE(i)*dB_vDE(j);
    end
end

N = zeros(nument,1);
temp_xi= zeros(nument,1);
temp_Beta= zeros(nument,1);

disp('element')
patchIndex
i
nument

PHTelem{patchIndex}(i).polyDegree

for t=1:nument
    
    if PHTelem{patchIndex}(i).polyDegree(t)==p
        N(t) = B*(PHTelem{patchIndex}(i).modifiedC(t,1:(p+1)*(q+1)))';
        temp_xi(t) = dBxi*(PHTelem{patchIndex}(i).modifiedC(t,1:(p+1)*(q+1)))';
        temp_Beta(t) =dBeta*(PHTelem{patchIndex}(i).modifiedC(t,1:(p+1)*(q+1)))';
    else
        N(t) = B_DE*(PHTelem{patchIndex}(i).modifiedC(t,:))';
        temp_xi(t) = dBxi_DE*(PHTelem{patchIndex}(i).modifiedC(t,:))';
        temp_Beta(t) = dBeta_DE*(PHTelem{patchIndex}(i).modifiedC(t,:))';
    end
    
end


%form B-spline basis functions using the Bezier extraction operator
% B = B*Ce';
% 
% dBxi = dBxi*Ce';
% dBeta = dBeta*Ce';
% 
% N = B;
% dN = [dBxi; dBeta];
dN = [temp_xi;temp_Beta];
