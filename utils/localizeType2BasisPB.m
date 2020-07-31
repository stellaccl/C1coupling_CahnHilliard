function [ PHUTelem] = localizeType2BasisPB(PHUTelem)
%remove the basis functions that have no support on a given element
%fix up and down nodes to remains periodic boundary


epsilon = 1e-10; %threshold for detecting empty support
list=[];
for indexPatch=1:length(PHUTelem)
    
    for indexElem=1:length(PHUTelem{indexPatch})
        
        if isempty(PHUTelem{indexPatch}(indexElem).children)
            
            modifiedC=PHUTelem{indexPatch}(indexElem).modifiedC;
            removeBasisIndex=[];
            for indexBasis=1:size(modifiedC,1)
                %nonzerosCoef=nonzeros(modifiedC(indexBasis,:));
                coefsNorm = norm(modifiedC(indexBasis,:));
                if coefsNorm<epsilon
                    removeBasisIndex=[ removeBasisIndex,indexBasis];
                end
            end
            
            if~isempty(removeBasisIndex)
            %    disp(['Patch ', num2str(indexPatch), ' element ', num2str(indexElem)])
            %    removeBasisIndex
                type2BasisIndex=PHUTelem{indexPatch}(indexElem).type2nodes;
                type2BasisIndex=sort(type2BasisIndex);
                remainType2BasisIndex=setdiff(type2BasisIndex,removeBasisIndex);
                %remainType2Basis=PHUTelem{indexPatch}(indexElem).nodesGlobal(remainType2BasisIndex);
                newNumType2Basis=length(remainType2BasisIndex);
                newType2BasisIndex=type2BasisIndex(1:newNumType2Basis);
                newRemoveBasis=type2BasisIndex(newNumType2Basis+1:end);
                
                tempGlobalBasis=PHUTelem{indexPatch}(indexElem).nodesGlobal;
                tempGlobalBasis(newType2BasisIndex)=tempGlobalBasis(remainType2BasisIndex);
                tempGlobalBasis(newRemoveBasis)=[];
                PHUTelem{indexPatch}(indexElem).nodesGlobal=tempGlobalBasis;
                
                % ====================================================================
                tempType2ModifiedC=PHUTelem{indexPatch}(indexElem).modifiedC;
                tempType2ModifiedC(newType2BasisIndex,:)=  tempType2ModifiedC(remainType2BasisIndex,:);
                tempType2ModifiedC(newRemoveBasis,:)=[];
                PHUTelem{indexPatch}(indexElem).modifiedC=tempType2ModifiedC;
                %==================================================================
 
            end
            
        end
    end
end

% tempGlobalBasis=PHUTelem{1}(8).nodesGlobal;
% saveBasis=tempGlobalBasis([3,4,7,8]);
% tempGlobalBasis([3,4,7,8])=tempGlobalBasis([11,12,15,16]);
% tempGlobalBasis([11,12,15,16])=saveBasis;
% PHUTelem{1}(8).nodesGlobal=tempGlobalBasis;
% 
% tempModifiedC=PHUTelem{1}(8).modifiedC;
% saveModifiedC=tempModifiedC([3,4,7,8],:);
% tempModifiedC([3,4,7,8],:)=tempModifiedC([11,12,15,16],:);
% tempModifiedC([11,12,15,16],:)=saveModifiedC;
% PHUTelem{1}(8).modifiedC=tempModifiedC;
% 
% tempGlobalBasis=PHUTelem{2}(7).nodesGlobal;
% saveBasis=tempGlobalBasis([1,2,5,6]);
% tempGlobalBasis([1,2,5,6])=tempGlobalBasis([9,10,13,14]);
% tempGlobalBasis([9,10,13,14])=saveBasis;
% PHUTelem{2}(7).nodesGlobal=tempGlobalBasis;
% 
% tempModifiedC=PHUTelem{2}(7).modifiedC;
% saveModifiedC=tempModifiedC([1,2,5,6],:);
% tempModifiedC([1,2,5,6],:)=tempModifiedC([9,10,13,14],:);
% tempModifiedC([9,10,13,14],:)=saveModifiedC;
% PHUTelem{2}(7).modifiedC=tempModifiedC;

end

