function [UniqueFams] = FindUniqueFamilies(cellOfVectors)
allVects = [];
for i = 1 : size(cellOfVectors,1)
    for j = i : size(cellOfVectors, 2)
        if ~isempty(cellOfVectors{i, j})
            numVectsij = length(cellOfVectors{i,j});
            for k = 1:numVectsij
                thisVect = cellOfVectors{i,j}{k};
                thisVect = sort(abs(thisVect));
                allVects = [allVects, thisVect];
            end
        end
    end
end

numAllVects = size(allVects, 2);
UniqueFams = {allVects(:, 1)};
for i = 1: numAllVects
    iVect = allVects(:, i);
    for j = (i+1): numAllVects
        isunique = 0;
        jVect = allVects(:, j);
        if all(abs(iVect - jVect) < 1e-5)
            for l = 1:length(UniqueFams)
                if ~(abs(UniqueFams{l} - jVect) < 1e-5)
                    isunique = 1;
                else
                    isunique = 0;
                    break;
                end
            end
        end
        if (isunique)
            UniqueFams = [UniqueFams; {jVect}];
        end    
    end
end

end