function[EI] = determine_EC_or_Ei(lambda, nu)


mu = 1.0;

lambda1 = 1 - lambda;

% solve for lambda2, the stretch in the x2 direction (use the positive
% solution).
lambda2_eq = [(mu + lambda1^2*2*nu*mu/(1 - nu)), -(2*mu*nu/(1 - nu))*lambda1,  -mu];
sol = roots(lambda2_eq);
lambda2Index = sol > 1.0;
lambda2 = sol(lambda2Index);

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
% L1212 = L2121;
L1221 = -2*lambda1*lambda2*dW_dII;
L2112 = L1221;

characteristic_eq = [L2121*L2222, 0,  (L2121^2 + L2222*L1111 - (L2112 + L2211)^2), 0, L1111*L2121];
phi_roots = roots(characteristic_eq);

if(all(abs(real(phi_roots)) < 1e-12))
    EI = 1;
else
    EI = 0;
end

end