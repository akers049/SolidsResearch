clear
clc
close all

stepSize = 0.015;

kappa_eval = 1.6;
w_vector = 0.1:stepSize:10;
solutionVector = zeros(size(w_vector));
doubleRootVector = zeros(size(w_vector));
initialGuessDoubleRoot = 0.53; %0.53; %0.5462; %0.53
initialGuess = 0.94;
for i = 1:length(w_vector)
    funHandle = @(x) evaluate_determinant_M(w_vector(i), kappa_eval, x);
    solutionVector(i) = fzero(funHandle, initialGuess);
    while (i > 2) && ...
            ((abs(2*solutionVector(i - 1) - solutionVector(i - 2) - solutionVector(i)) > stepSize/2) ...
            || (solutionVector(i) > 1.0))
        initialGuess = initialGuess + stepSize/10*(rand - 0.5);
        solutionVector(i) = fzero(funHandle, initialGuess);
    end
    
    doubleRootVector(i) = fzero(funHandle, initialGuessDoubleRoot);
    while (i > 2) && ...
            ((abs(2*doubleRootVector(i - 1) - doubleRootVector(i - 2) - doubleRootVector(i)) > stepSize/4) ...
            || (doubleRootVector(i) > 1.0))
        initialGuessDoubleRoot = initialGuessDoubleRoot + stepSize/10*(rand - 0.5);
        doubleRootVector(i) = fzero(funHandle, initialGuessDoubleRoot);
    end
    if (i > 1)
        initialGuess =  2*solutionVector(i) - solutionVector(i - 1);
        initialGuessDoubleRoot = (doubleRootVector(i) - doubleRootVector(i -1)) + doubleRootVector(i);
    else
        initialGuess = solutionVector(i) - stepSize*10;
        initialGuessDoubleRoot = doubleRootVector(i);
    end
end

det_M_tilde = zeros(size(w_vector));
for i = 1:length(w_vector)
    det_M_tilde(i) = evaluate_determinant_Mtilde(w_vector(i) , kappa_eval, doubleRootVector(i));
   if abs(det_M_tilde(i)) < stepSize*10
       solutionVector(i) = doubleRootVector(i); 
   end
end

plot(w_vector, solutionVector, 'LineWidth', 1.3)
hold on;
plot(w_vector, doubleRootVector, 'LineWidth', 1.3)
grid on
ylabel('\lambda_{c}      ', 'rot', 1)
xlabel('\omega', 'fontweight', 'bold', 'fontsize', 13)
title('Ciritcal Load for Compressed Strip')
axis([0 8 0, 1]);
