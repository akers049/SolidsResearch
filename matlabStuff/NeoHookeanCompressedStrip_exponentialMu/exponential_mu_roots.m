function[systemRoots] = exponential_mu_roots(lambda, nu, kappa, w)
mu0 = 1.0;


lambda1 = 1 - lambda;

% solve for lambda2, the stretch in the x2 direction (use the positive
% solution).
lambda2_eq = [(mu0 + lambda1^2*2*nu*mu0/(1 - nu)), -(2*mu0*nu/(1 - nu))*lambda1,  -mu0];
sol = roots(lambda2_eq);
lambda2Index = sol > 1.0;
lambda2 = sol(lambda2Index);

% calculalate derivatives of energy with the invariants I and II (trace and
% determinant)

dW_dI = mu0/2;
dW_dII = -mu0/(2*lambda1^2*lambda2^2) + (mu0*nu/(1- nu))*(1 - 1/lambda1/lambda2);
d2W_dIdI = 0;
d2W_dIdII = 0;
d2W_dIIdII = mu0/(2*lambda1^4*lambda2^4) + ...
    (mu0*nu/(2*(1-nu)))*(1/(lambda1^3*lambda2^3));

% calculate components of the incremental moduli along principal branch
L1111 = 4*lambda1^2*(d2W_dIdI + 2*lambda2^2*d2W_dIdII + lambda2^4*d2W_dIIdII) + ...
    2*(dW_dI + lambda2^2*dW_dII);
L2222 = 4*lambda2^2*(d2W_dIdI + 2*lambda1^2*d2W_dIdII + lambda1^4*d2W_dIIdII) + ...
    2*(dW_dI + lambda1^2*dW_dII);
L1122 = 4*lambda1*lambda2*(d2W_dIdI + (lambda1^2 + lambda2^2)*d2W_dIdII + ...
    lambda1^2*lambda2^2*d2W_dIIdII + dW_dII);
L2121 = 2*dW_dI;
L1212 = L2121;
L1221 = -2*lambda1*lambda2*dW_dII;
L2112 = L1221;

a = L2222*L1212;
b = 2*kappa*L2222*L1212;
c = w^2*(L1122^2 + 2*L1122*L2112 - L1212^2 + L2222*L1212*(kappa^2)/(w^2) + L2112^2 - L1111*L2222);
d = kappa*w^2*(L1122^2 + 2*L2112*L1122 - L1212^2 + L2112^2 - L1111*L2222);
e = w^2*(L1122*L2112*kappa^2 + L1111*L1212*w^2);

characterist_eq = [a b c d e];
systemRoots = roots(characterist_eq);

end