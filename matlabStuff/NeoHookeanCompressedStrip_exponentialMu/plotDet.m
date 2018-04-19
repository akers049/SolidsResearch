
kappa_eval = 1.0;
w_eval = 1;
lambdas = 0.36:0.0001:0.37;
detVals = zeros(size(lambdas));
% syms lambda;
% symFunction = symfun(m(w_eval, kappa_eval, lambda), lambda);

%%
for i = 1:length(lambdas)
    detVals(i) = evaluate_determinant_fast(w_eval, kappa_eval, lambdas(i));
end

plot(lambdas, detVals)
