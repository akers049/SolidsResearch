function [mHats, bs] = SolveTwinningAndInterfaceEqs(Us)

% We will store the a vectors for the corresponding I and J variants of
% martensite into aVectors(I,J). aVectors should be symetric. (I think)
% Same for the nHat vectors.

% get the number of variants
numVariants = length(Us);

% allocate the cells to store the data
aVectors = cell(numVariants, numVariants);
nHatVectors = cell(numVariants, numVariants);

% Loop over all the unlike pairs of variants
for i = 1:numVariants
    for j = (i+1):numVariants
        % Pull the transformation matricies for variants I and J from the
        % cell array of transformation matricies.
        UI = Us{i};
        UJ = Us{j};
        
        % Calculate the C matrix
        C = inv(transpose(UJ))*transpose(UI)*UI*inv(UJ);
        
        % Solution only if the C matrix is not the identity matrix
        if ~(isequal((abs(C - eye(3)) < 1e-6), ones(3,3)))
            % find and order the eigenvalues of C. 
            [lamdas, eHats] = calcAndOrderEigs(C);

            % check that the eigenvalues stradle 1, and that lamda(2) is equal
            % to 1.
            if  ((abs(lamdas(2) - 1) < 1e-6) && (lamdas(1) <= 1) && (lamdas(3) >= 1));
                % Now we can calculate the two solutions, corresponding to if
                % kappa=-1 or kappa=1. We will calculate these in the function
                % 'calcTwinVectors'
                [a1, nHat1] = calcTwinVectors(lamdas, eHats, UJ, 1);
                [a2, nHat2] = calcTwinVectors(lamdas, eHats, UJ, -1);
                
                % Now we just store the data in the cell arrays
                aVectors{i,j} = {a1; a2};
                aVectors(j,i) = aVectors(i,j);

                nHatVectors{i,j} = {nHat1; nHat2};
                nHatVectors{j,i} = nHatVectors{i,j};
            end
        end
    end
end
%%
% Now that we have the solutions to the twining equation, we will now
% calculate all of the corresponding solutions to the austensite-martensite
% interface equation. We will store the b vectors into cell array bs, and 
% mHat vectors into a cell array, mHatVectors. mHatVectors(I,J) contains
% the mHat vectors that correspond to the twinning between variants I and J
% of martensite.

% Alocate cell arrays to store the mHat and b vectors
mHats = cell(numVariants, numVariants);
bs = cell(numVariants, numVariants);

% loop over unlike pairs of variants
for i = 1:numVariants;
    for j = (i+1):numVariants;
        % indx is used for soring the possibly variable number of interface
        % solutions
        indx = 1;
        
        % UI is actually the J variant used for the twinning equation...
        % Tricky! (equation 7.9 page 112)
        UI = Us{j};
        
        % Loop over the (possibly) two solutions to the twinning equation
        for k = 1:2
           % if there were no solutions to the twinning equation, then just
           % break and move to next pair of variants
           if (isempty(aVectors{i, j}))
               break;
           end
           % but if there are solutions, there will be two, and we will
           % loop over both of them
           
           % pull the solutions for the twinning equation from the cell
           % arrays (a and nHat)
           a = aVectors{i,j}{k};
           nHat = nHatVectors{i,j}{k};
           
           % Calculate the delta and nu (step 1 page 113, eq 7.13)
           [delta, nu] = calcDeltaNu(a, nHat, UI);
           
           % If the delta and nu meet the criteria of eq 7.14, then
           % continue
           if (delta <= -2 && nu >= 0)
               % calculate lamda
               lam = 0.5*(1 - sqrt(1 + 2/delta));
               
               % execute step 3 on page 113
               [b, mHat] = Step3Page113(lam ,a, nHat, UI);
               
               % store the mHats and bs. Use indx because there could be a
               % variable number of solutions.
               mHats{i, j}(indx:(indx+1), 1) = mHat;
               bs{i, j}(indx:(indx+1), 1) = b;
               indx = indx+2;
               
               % Now if delta is less than 2, we have to replace lamda with
               % (1 - lamda), and do step 3 again.
               if ((delta + 1e-6 < 2))
                   lam = (1- lam);
                   [b, mHat] = Step3Page113(lam ,a, nHat, UI);
                   mHats{i, j}(indx:(indx+1), 1) = mHat;
                   bs{i, j}(indx:(indx+1), 1) = b;
                   indx = indx+2;
               end
           end
        end
    end
end

end
        
        