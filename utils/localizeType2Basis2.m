function [ PHUTelem] = localizeType2Basis2(PHUTelem,type2Basis)
disp('localizing type 2 basis function')
%remove the basis functions that have no support on a given element
%relocate according to type 2 basis function location

epsilon = 1e-10; %threshold for detecting empty support

for indexPatch=1:length(PHUTelem)
    for indexElem=1:length(PHUTelem{indexPatch})
        
        if isempty(PHUTelem{indexPatch}(indexElem).children)
            modifiedC=PHUTelem{indexPatch}(indexElem).modifiedC;
            removeBasis=[];
            for indexBasis=1:size(modifiedC,1)
                
                if ismember(indexBasis,type2Basis)  
                    %nonzerosCoef=nonzeros(modifiedC(indexBasis,:));
                    coefsNorm = norm(modifiedC(indexBasis,:));
                    if coefsNorm<epsilon
                        removeBasis=[ removeBasis,indexBasis];
                    end
                    
                end
                
            end
            
            if~isempty(removeBasis)
                %  disp(['Patch ', num2str(indexPatch), ' element ', num2str(indexElem)])
                %  removeBasis
                tempGlobalBasis=PHUTelem{indexPatch}(indexElem).nodesGlobal;
                tempGlobalBasis(removeBasis)=[];
                
                PHUTelem{indexPatch}(indexElem).nodesGlobal=tempGlobalBasis;
                
                tempModifiedC=PHUTelem{indexPatch}(indexElem).modifiedC;
                tempModifiedC(removeBasis,:)=[];
                PHUTelem{indexPatch}(indexElem).modifiedC=tempModifiedC;
                % disp([num2str(length(PHUTelem{indexPatch}(indexElem).nodesGlobal)), ' basis remaining'])
            end
            
        end
    end
end



end

