function [criticalLoad, criticalWavelength] = calculateFirstBifLoad_pieceConst(L1, kappa, initGuess)

kappa_eval = kappa;
stepSize = 0.1;
wavelength_vector = 0.1:stepSize:60;
solutionVector = zeros(size(wavelength_vector));
initialGuess = initGuess;
onSolution = 1 ;
errorFlag = 0;
for i = 1:length(wavelength_vector)
    funHandle = @(x) piece_const_mu_evaluate_determinant_M(L1, kappa_eval, wavelength_vector(i), x);

    solutionVector(i) = fzero(funHandle, initialGuess);

    count = 0;
    guessIncrement = 0.0;
    if (onSolution == 0)
        while (abs(solutionVector(i - 1) - solutionVector(i)) > stepSize)
            count = count + 1;
            if (mod(count, 2) == 0)
              guessIncrement = guessIncrement + stepSize/10; %*(rand - 0.5);
              solutionVector(i) = fzero(funHandle, (initialGuess + guessIncrement));
            else
              solutionVector(i) = fzero(funHandle, abs(initialGuess - guessIncrement));
            end

            if (count > 10/stepSize)
                fprintf('%d super fucked\n', i);
                errorFlag = 1;
                break;
            end
        end
        onSolution = 1;
    else
        while (i > 2) && ...
                ((abs(2*solutionVector(i - 1) - solutionVector(i - 2) - solutionVector(i)) > stepSize/20) ...
                || (solutionVector(i) > 1.0)) && (onSolution == 1)
            count = count + 1;
            if (mod(count, 2) == 0)
              guessIncrement = guessIncrement + stepSize/20; %*(rand - 0.5);
              solutionVector(i) = fzero(funHandle, (initialGuess + guessIncrement));
            else
              solutionVector(i) = fzero(funHandle, abs(initialGuess - guessIncrement));
            end

            if (count > 20/stepSize)
                solutionVector(i) = initialGuess;
                fprintf('%d really fucked\n', i);
                onSolution = 0;
                break;
            end
        end
    end
    if (errorFlag == 1)
        break
    end
    
    if (i > 1)
        initialGuess =  2*solutionVector(i) - solutionVector(i - 1);
    else
        initialGuess = solutionVector(i);
    end
end


plot(wavelength_vector, solutionVector, 'LineWidth', 1.3)
grid on
ylabel('\lambda_{c}      ', 'rot', 1)
xlabel('Wavelength', 'fontweight', 'bold', 'fontsize', 13)
title('Ciritcal Load for Compressed Strip')

criticalLoad = min(solutionVector);
criticalWavelength =  wavelength_vector(find(min(solutionVector) == solutionVector));

end