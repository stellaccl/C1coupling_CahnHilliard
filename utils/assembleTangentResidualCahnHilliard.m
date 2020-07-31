function [ tangent, residual ] = assembleTangentResidualCahnHilliard( PHTelem, basisFun, sizeBasis, p, q,C,Cdot,theta,alpha,shift)
%assembles tangent and residual for CahnHilliard equation
%uses saved basis functions

%Gauss points
ngauss_x = p+1;
ngauss_y = q+1;
[gauss_weight_x] = quadrature( ngauss_x, 'GAUSS', 1 );
[gauss_weight_y] = quadrature( ngauss_y, 'GAUSS', 1 );

%initialize residual vector
residual= zeros(sizeBasis,1);

%initialize tangent matrix
%allocate memory for the triplet arrays
indexCounter = 0;
for patchIndex=1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            nument = size(PHTelem{patchIndex}(i).C,1);
            indexCounter = indexCounter + nument^2;
        end
    end
end
%indexCounter
II = zeros(1,indexCounter);
JJ = zeros(1,indexCounter);
S = zeros(1, indexCounter);

indexCounter = 0;

%%tangent = sparse(sizeBasis,sizeBasis);

for patchIndex = 1:length(PHTelem)
    %     patchIndex
    %     pause
    for elemIndex=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(elemIndex).children)
            xmin = PHTelem{patchIndex}(elemIndex).vertex(1);
            xmax = PHTelem{patchIndex}(elemIndex).vertex(3);
            ymin = PHTelem{patchIndex}(elemIndex).vertex(2);
            ymax = PHTelem{patchIndex}(elemIndex).vertex(4);
            
            %the jacobian of the transformation from [-1,1]x[-1,1] to
            %[xmin, xmax]x [ymin, ymax]
            scalefac = (xmax - xmin)*(ymax - ymin)/4;
            nument = size(PHTelem{patchIndex}(elemIndex).modifiedC,1);
            sctrx = PHTelem{patchIndex}(elemIndex).nodesGlobal(1:nument);
            
            localTangent = zeros(nument, nument);
            localResidual = zeros(nument, 1);
            
            tempC=C(PHTelem{patchIndex}(elemIndex).nodesGlobal);
            tempCdot=Cdot(PHTelem{patchIndex}(elemIndex).nodesGlobal);
            %elemIndex
            %loop over the ngauss_x x ngauss_y gauss points on each element
            for jj=1:ngauss_y
                for ii=1:ngauss_x
                    %evaluate the derivatives of the mapping from parameter
                    %space to physical space
                    
                    R = basisFun.R{patchIndex,elemIndex,ii,jj};
                    dRdx = basisFun.dRdx{patchIndex,elemIndex,ii,jj};
                    dRdy = basisFun.dRdy{patchIndex,elemIndex,ii,jj};
                    d2Rdx = basisFun.d2Rdx{patchIndex,elemIndex,ii,jj};
                    d2Rdy = basisFun.d2Rdy{patchIndex,elemIndex,ii,jj};
                    d2Rdxdy = basisFun.d2Rdxdy{patchIndex,elemIndex,ii,jj};
                    
                    dxdxi = basisFun.dxdxi{patchIndex,elemIndex,ii,jj};
                    d2xdxi2 = basisFun.d2xdxi2{patchIndex,elemIndex,ii,jj};
                    dR = basisFun.dR{patchIndex,elemIndex,ii,jj};
                    
                    laplacian = basisFun.laplacian{patchIndex,elemIndex,ii,jj};
                    J1    = basisFun.J1{patchIndex,elemIndex,ii,jj};
                    
                    %tangent part
                    [localC, dLocalC, del2_c, localC_t]=computeConcentrationInfo2(tempC,tempCdot,R,dRdx,dRdy,d2Rdx,d2Rdy,d2Rdxdy,dxdxi,d2xdxi2);
                    [~,dmu,d2mu ] = chemicalPotential( localC, theta,alpha );
                    [M,dM,d2M] = mobility( localC);
                    t1=M*dmu+dM*del2_c;
                    t2=dM*dmu+M*d2mu+d2M*del2_c;
                                                            
                    localTangentTemp = computeLocalTanCH2(R,dR,laplacian,shift,dM,t1,t2,dLocalC,del2_c,M);                    
                    localTangent = localTangent + localTangentTemp*scalefac*gauss_weight_x(ii).*gauss_weight_y(jj).*J1;
                    
                    %residual part
                    localResidualTemp=R*localC_t+( dLocalC*dR)'*t1+(laplacian)*M*del2_c;
                    
                    localResidualTemp=localResidualTemp.*scalefac.*gauss_weight_x(ii).*gauss_weight_y(jj).*J1;
                    localResidual = localResidual + localResidualTemp;
                    
                end
            end
            
            %localTangent = localTangent+localTangent'-diag(diag(localTangent));
            residual(sctrx) =  residual(sctrx) + localResidual;
            
            %localTangent
            %tangent(scrtx,scrtx) = tangent(scrtx,scrtx) + localTangent;
            II(indexCounter+1:indexCounter+nument^2) = repmat(sctrx,1,nument);
            JJ(indexCounter+1:indexCounter+nument^2) = reshape(repmat(sctrx',1,nument)',1,nument^2);
            S(indexCounter+1:indexCounter+nument^2) = reshape(localTangent,1,nument^2);
            indexCounter = indexCounter + nument^2;
        end
    end
end


tangent = sparse(II, JJ, S, sizeBasis, sizeBasis);

end

