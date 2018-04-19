function [lamdas, eHats] = calcAndOrderEigs(C)
% This function takes a matrix C and finds its eigenvalues and
% eigenvectors, which it puts in order. eigen values are ordered from least
% to greatest going from top to bottom of a column vector lamdas. Their
% corresponding eigenvectors are the columns from left to right of a square
% matrix eHats.

[eigVects , eigVals] = eig(C);
eigVals = eigVals*ones(size(C,1), 1);
lamdas = sort(eigVals);
eHats = zeros(size(C));
for k = 1:size(C,1);
    eHats(:, k) = eigVects(:, find(lamdas(k) == eigVals));
end
        
end