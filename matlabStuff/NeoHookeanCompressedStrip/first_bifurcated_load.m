clear
clc
close all

mu = 1;
nu = 0.33;
L = 1;

% define bifucation variable lambda, and a dummy variable x
syms lambda x;

% give lambda1, the stretch in the x1 direction
lambda1(lambda) = 1 - lambda;

% solve for lambda2, the stretch in the x2 direction (use the positive
% solution).
lambda2_eq = (mu + lambda1^2*2*nu*mu/(1 - nu))*x^2 - (2*mu*nu/(1 - nu))*lambda1*x - mu == 0;
sol = solve(lambda2_eq, x);
solFun(lambda) = sol;
lambda2Index = find(double(solFun(0.3)) > 0);
lambda2(lambda) = (sol(lambda2Index));

% calculalate derivatives of energy with the invariants I and II (trace and
% determinant)

dW_dI = mu/2;
dW_dII = -mu/(2*lambda1^2*lambda2^2) + (mu*nu/(1- nu))*(1 - 1/lambda1/lambda2);
d2W_dIdI = 0;
d2W_dIdII = 0;
d2W_dIIdII = mu/(2*lambda1^4*lambda2^4) + ...
    (mu*nu/(2*(1-nu)))*(1/(lambda1^3*lambda2^3));

% calculate components of the incremental moduli along principal branch
L1111 = 4*lambda1^2*(d2W_dIdI + 2*lambda2^2*d2W_dIdII + lambda2^4*d2W_dIIdII) + ...
    2*(dW_dI + lambda2^2*dW_dII);
L2222 = 4*lambda2^2*(d2W_dIdI + 2*lambda1^2*d2W_dIdII + lambda1^4*d2W_dIIdII) + ...
    2*(dW_dI + lambda1^2*dW_dII);
L1122 = 4*lambda1*lambda2*(d2W_dIdI + (lambda1^2 + lambda2^2)*d2W_dIdII + ...
    lambda1^2*lambda2^2*d2W_dIIdII + dW_dII);
L2211 = L1122;
L2121 = 2*dW_dI;
L1212 = L2121;
L1221 = -2*lambda1*lambda2*dW_dII;
L2112 = L1221;

characteristic_eq = (L2121*L2222*x^2 + (L2121^2 + L2222*L1111 - (L2112 + L2211)^2)*x + L1111*L2121) == 0;
phi_squared = solve(characteristic_eq, x);

alpha(lambda) = simplify(sqrt(-phi_squared(1)));
beta(lambda) = simplify(sqrt(-phi_squared(2)));

K_excee(x) = (L1111 - L1212*x^2)/((L2112 + L2211)*x);
K_alpha(lambda) = simplify(K_excee(alpha));
K_beta(lambda) = simplify(K_excee(beta));

% now lets put in the equation to solve for
syms w;
toSolve(w, lambda) = tanh(alpha*w*L)*(L2112*K_alpha + L1212*alpha)*(L2222*beta*K_beta - L1122) - ...
     tanh(beta*w*L)*(L2112*K_beta + L1212*beta)*(L2222*alpha*K_alpha - L1122);
 
%% 
% lambdas = 0.0001:0.001:0.97;
% tmp_vect = zeros(size(lambda));
% for i = 1: length(lambdas)
%    tmp_vect(i) = double(toSolve(0.11, lambdas(i)));
% end
% plot(lambdas, tmp_vect);

w_vector = 0.11:0.1:8;
solutionVector = zeros(size(w_vector));
initialGuess = 0.9581;
for i = 1:length(w_vector)
    solutionVector(i) = vpasolve((toSolve(w_vector(i), lambda) ==0), lambda, initialGuess);
    initalGuess = solutionVector(i);
end

plot(w_vector, solutionVector, 'LineWidth', 1.3)
grid on
ylabel('\lambda_{c}      ', 'rot', 1)
xlabel('\omega', 'fontweight', 'bold', 'fontsize', 13)
title('Ciritcal Load for Compressed Strip')






