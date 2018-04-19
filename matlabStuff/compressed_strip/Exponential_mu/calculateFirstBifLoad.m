function [criticalLoad, criticalWavelength] = calculateFirstBifLoad(kappa, initGuess)

kappa_eval = kappa;
stepSize = 0.1;
wavelength_vector = 0.3:stepSize:30;
solutionVector = zeros(size(wavelength_vector));
initialGuess = initGuess;
onSolution = 1 ;
fuckedFlag = 0.0;
for i = 1:length(wavelength_vector)
    funHandle = @(x) evaluate_determinant_M(wavelength_vector(i), kappa_eval, x);
    solutionVector(i) = fzero(funHandle, initialGuess);
    
    count = 0;
    guessIncrement = 0.0;
    if (onSolution == 0)
        while (abs(solutionVector(i - 1) - solutionVector(i)) > stepSize/1.0)
            count = count + 1;
            if (mod(count, 2) == 0)
              guessIncrement = guessIncrement + stepSize/100; %*(rand - 0.5);
              solutionVector(i) = fzero(funHandle, (initialGuess + guessIncrement));
            else
              solutionVector(i) = fzero(funHandle, abs(initialGuess - guessIncrement));
            end

            if (count > 100/stepSize)
                fprintf('%d super fucked\n', i);
                fuckedFlag = 1.0;
                break;
            end
        end
        onSolution = 1;
    else
        while (i > 2) && ...
                ((abs(2*solutionVector(i - 1) - solutionVector(i - 2) - solutionVector(i)) > stepSize/100) ...
                || ((solutionVector(i) > 1.0) || (solutionVector(i) < 0.0)))
            count = count + 1;
            if (mod(count, 2) == 0)
              guessIncrement = guessIncrement + stepSize/100; %*(rand - 0.5);
              solutionVector(i) = fzero(funHandle, (initialGuess + guessIncrement));
            else
              solutionVector(i) = fzero(funHandle, abs(initialGuess - guessIncrement));
            end

            if (count > 100/stepSize)
                solutionVector(i) = initialGuess;
                fprintf('%d really fucked\n', i);
                onSolution = 0;
                break;
            end
        end
    end
    
    if (i > 1)
        initialGuess =  2*solutionVector(i) - solutionVector(i - 1);
    else
        initialGuess = solutionVector(i);
    end
    
    if (fuckedFlag == 1.0)
        break;
    end
end


plot(wavelength_vector, solutionVector, 'LineWidth', 1.3)
grid on
ylabel('\lambda_{c}      ', 'rot', 1)
xlabel('Wavelength', 'fontweight', 'bold', 'fontsize', 13)
titleStr = sprintf(['Ciritcal Load for exponential $\\mu$ \n $\\kappa$ = ', num2str(kappa)]);
title(titleStr, 'Interpreter', 'latex');
criticalLoad = min(solutionVector);
criticalWavelength =  wavelength_vector(find(min(solutionVector) == solutionVector));

end