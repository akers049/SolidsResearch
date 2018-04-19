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
alphas = root(characteristic_eq == 0, z);

% now can solve for B's:
B = sym(size(alphas));
for i = 1:4
    B(i) = (L1111*w^2 - kappa*L1212*alphas(i) - L1212*alphas(i)^2)/(w*alphas(i)*(L2112 + L1122) + w*kappa*L2112);
   % B(i) = (w*alphas(i)*(L1221 + L1122) + kappa*L1122*w)/(kappa*L2222*alphas(i) + L2222*alphas(i)^2 - L2121*w^2);
end

% now construct the system matrix:
syms systemMat;

systemMat = sym(zeros(4,4));
% first row
for i = 1:4
    systemMat(1, i) = B(i);
end

% second row
for i = 1:4
    systemMat(2, i) = L1212*alphas(i) + w*L2112*B(i);
end

% third row
for i = 1:4
    systemMat(3, i) = L2222*alphas(i)*B(i)*exp(1)^(alphas(i)*L) - w*L1122*exp(1)^(alphas(i)*L);
end

% fourth row
for i = 1:4
    systemMat(4, i) = w*L2112*B(i)*exp(1)^(alphas(i)*L) + L1212*alphas(i)*exp(1)^(alphas(i)*L);
end


m(w, kappa, lambda) = systemMat;


% toSolve(w, lambda) = m(1,1)*m(2,2)*m(3,3)*m(4,4) - m(1,1)*m(2,2)*m(3,4)*m(4,3) ...
%    - m(1,1)*m(2,3)*m(3,2)*m(4,4) + m(1,1)*m(2,3)*m(3,4)*m(4,2) + m(1,1)*m(2,4)*m(3,2)*m(4,3) ...
%    - m(1,1)*m(2,4)*m(3,3)*m(4,2) - m(1,2)*m(2,1)*m(3,3)*m(4,4) + m(1,2)*m(2,1)*m(3,4)*m(4,3) ...
%    
% + m(1,2)*m(2,3)*m(3,1)*m(4,4) - m(1,2)*m(2,3)*m(3,4)*m(4,1) - m(1,2)*m(2,4)*m(3,1)*m(4,3) ...
%    + m(1,2)*m(2,4)*m(3,3)*m(4,1) + m(1,3)*m(2,1)*m(3,2)*m(4,4) - m(1,3)*m(2,1)*m(3,4)*m(4,2) ...
%    - m(1,3)*m(2,2)*m(3,1)*m(4,4) + m(1,3)*m(2,2)*m(3,4)*m(4,1) + m(1,3)*m(2,4)*m(3,1)*m(4,2) ...
%    - m(1,3)*m(2,4)*m(3,2)*m(4,1) - m(1,4)*m(2,1)*m(3,2)*m(4,3) + m(1,4)*m(2,1)*m(3,3)*m(4,2) ...
%    + m(1,4)*m(2,2)*m(3,1)*m(4,3) - m(1,4)*m(2,2)*m(3,3)*m(4,1) - m(1,4)*m(2,3)*m(3,1)*m(4,2) ...
%    + m(1,4)*m(2,3)*m(3,2)*m(4,1);

%% Solving!!!!
kappa_eval = 1.3;

w_vector = 0.1:0.1:8;
solutionVector = zeros(size(w_vector));
initialGuess = 0.95;
for i = 1:length(w_vector)
    symMat = symfun(m(w_vector(i), kappa_eval, lambda), lambda);
    tic;
    solutionVector(i) = bisect_symbolic_determinant_root(symMat, (initialGuess - 0.1) , 0.96, 0.0001, 20);
    toc;
    clear symFunction;
    initialGuess = solutionVector(i);
end

plot(w_vector, solutionVector, 'LineWidth', 1.3)
grid on
ylabel('\lambda_{c}      ', 'rot', 1)
xlabel('\omega', 'fontweight', 'bold', 'fontsize', 13)
title('Ciritcal Load for Compressed Strip')


