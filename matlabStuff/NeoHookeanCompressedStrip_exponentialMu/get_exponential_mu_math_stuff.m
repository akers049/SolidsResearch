clear
clc
close all

mu0 = 1;
nu = 0.33;
L = 1;

% define bifucation variable lambda, and a dummy variable x
syms lambda x kappa;

% give lambda1, the stretch in the x1 direction
lambda1(lambda) = 1 - lambda;

% solve for lambda2, the stretch in the x2 direction (use the positive
% solution).
lambda2_eq = (mu0 + lambda1^2*2*nu*mu0/(1 - nu))*x^2 - (2*mu0*nu/(1 - nu))*lambda1*x - mu0 == 0;
sol = solve(lambda2_eq, x);
solFun(lambda) = sol;
lambda2Index = find(double(solFun(0.3)) > 0);
lambda2(lambda) = (sol(lambda2Index));

% calculalate derivatives of energy with the invariants I and II (trace and
% determinant)

dW_dI = mu0/2;
dW_dII = -mu0/(2*lambda1^2*lambda2^2) + (mu0*nu/(1- nu))*(1 - 1/lambda1/lambda2);
d2W_dIdI = 0;
d2W_dIdII = 0;
d2W_dIIdII = mu0/(2*lambda1^4*lambda2^4) + ...
    (mu0*nu/(2*(1-nu)))*(1/(lambda1^3*lambda2^3));

% calculate components of the incremental moduli along principal branch
L1111(lambda) = 4*lambda1^2*(d2W_dIdI + 2*lambda2^2*d2W_dIdII + lambda2^4*d2W_dIIdII) + ...
    2*(dW_dI + lambda2^2*dW_dII);
L2222(lambda) = 4*lambda2^2*(d2W_dIdI + 2*lambda1^2*d2W_dIdII + lambda1^4*d2W_dIIdII) + ...
    2*(dW_dI + lambda1^2*dW_dII);
L1122(lambda) = 4*lambda1*lambda2*(d2W_dIdI + (lambda1^2 + lambda2^2)*d2W_dIdII + ...
    lambda1^2*lambda2^2*d2W_dIIdII + dW_dII);
L2211(lambda) = L1122;
L2121 = 2*dW_dI;
L1212 = L2121;
L1221(lambda) = -2*lambda1*lambda2*dW_dII;
L2112(lambda) = L1221;

% Make the characteristic equation;
syms a b c d e w z;

a = simplify(L2222*L1212);
b = simplify(2*kappa*L2222*L1212);
c = simplify(w^2*(L1122^2 + 2*L1122*L2112 + L2112^2 - L1212^2  - L1111*L2222) + L2222*L1212*(kappa^2));
d = simplify(kappa*w^2*(L1122^2 + 2*L2112*L1122  + L2112^2 - L1212^2 - L1111*L2222));
e = simplify(w^2*(L1122*L2112*kappa^2 + L1111*L1212*w^2));

characteristic_eq = a*z^4 + b*z^3 + c*z^2 + d*z + e;

syms alphas B;
% solve for alphas
alphas_sym(w, kappa, lambda) = root(characteristic_eq == 0, z);

save('exponential_mu_math_stuff.mat', 'L', 'lambda', 'kappa', 'w', 'alphas_sym', 'L1111', 'L2222', 'L1122', 'L1212', 'L2112'); 


