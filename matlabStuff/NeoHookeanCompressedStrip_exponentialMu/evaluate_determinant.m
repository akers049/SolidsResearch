function [detVal] = evaluate_determinant(w_eval, kappa_eval, lambda_eval)
load('exponential_mu_math_stuff.mat');

alphas_double = double(alphas_sym(w_eval, kappa_eval, lambda_eval));

B = zeros(4, 1);
for i = 1:4
    B(i) = (double(L1111(lambda_eval))*w_eval^2 - ...
        kappa_eval*L1212*alphas_double(i) - ...
        L1212*alphas_double(i)^2) / ...
        (w_eval*alphas_double(i)*double(L2112(lambda_eval) + ...
        L1122(lambda_eval)) + w_eval*kappa_eval*double(L2112(lambda_eval)));
end

systemMat = zeros(4,4);
for i = 1:4
    systemMat(1, i) = B(i);
end

% second row
for i = 1:4
    systemMat(2, i) = L1212*alphas_double(i) + ...
        w_eval*double(L2112(lambda_eval))*B(i);
end

% third row
for i = 1:4
    systemMat(3, i) = double(L2222(lambda_eval))*alphas_double(i)*B(i)*exp(alphas_double(i)*L) - ...
        w_eval*double(L1122(lambda_eval))*exp(alphas_double(i)*L);
end

% fourth row
for i = 1:4
    systemMat(4, i) = w_eval*double(L2112(lambda_eval))*B(i)*exp(alphas_double(i)*L) + ...
        L1212*alphas_double(i)*exp(alphas_double(i)*L);
end

detVal = real(det(systemMat));
end