
kappa_eval = 3.4;
w_eval = 2*pi/1.4; %2*pi/1.7;
lambdas = 0.1:0.0001:0.7;
detVals = zeros(1, length(lambdas));
% syms lambda;
% symFunction = symfun(m(w_eval, kappa_eval, lambda), lambda);

%%
for i = 1:length(lambdas)
    detVals(:, i) = evaluate_determinant_M(w_eval, kappa_eval, lambdas(i));
end
figure(2)
hold on
for i = 1:size(detVals, 1)
   plot(lambdas, detVals(i, :), '.', 'Markersize', 5) 
end
% axis([lambdas(1), lambdas(end), -1 1])
