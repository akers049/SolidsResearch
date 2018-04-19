clear
clc
close all

kappa_eval = 1.0;
w_vector = 0.1:0.01:10;
solutionVector = zeros(size(w_vector));
initialGuess = 0.94;
for i = 1:length(w_vector)
    funHandle = @(x) evaluate_determinant_fast(w_vector(i), kappa_eval, x);
    solutionVector(i) = fzero(funHandle, initialGuess);
    while (i > 1) &&(abs(solutionVector(i) - solutionVector(i - 1)) > 0.01)
        initialGuess = initialGuess + 0.01*(rand - 0.5);
        solutionVector(i) = fzero(funHandle, initialGuess);
    end
    initialGuess = solutionVector(i) + 0.001;
end

plot(w_vector, solutionVector, 'LineWidth', 1.3)
grid on
ylabel('\lambda_{c}      ', 'rot', 1)
xlabel('\omega', 'fontweight', 'bold', 'fontsize', 13)
title('Ciritcal Load for Compressed Strip')
axis([0 8 0.3, 1]);
