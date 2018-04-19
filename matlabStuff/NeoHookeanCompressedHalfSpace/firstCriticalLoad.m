function [ lambda_c, wavenumber ] = firstCriticalLoad( mu, nu )

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
lambda2(lambda) = sol(lambda2Index);

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
%L1111 = simplify(L1111);
L2222 = 4*lambda2^2*(d2W_dIdI + 2*lambda1^2*d2W_dIdII + lambda1^4*d2W_dIIdII) + ...
    2*(dW_dI + lambda1^2*dW_dII);
%L2222 = simplify(L2222);
L1122 = 4*lambda1*lambda2*(d2W_dIdI + (lambda1^2 + lambda2^2)*d2W_dIdII + ...
    lambda1^2*lambda2^2*d2W_dIIdII + dW_dII);
%L1122 = simplify(L1122);
L2211 = L1122;
L2121 = 2*dW_dI;
L1212 = L2121;
L1221 = -2*lambda1*lambda2*dW_dII;
%L1221 = simplify(L1221);
L2112 = L1221;


% solve the biquadratic for z1 and z2 (solutions with negative imaginary
% parts)
% syms z;
% biQuadratic_zEquation = (L1111 + z^2*L1212)*(L2121+z^2*L2222) - z^2*(L1122 + L1221)^2 == 0;
% solz = solve(biQuadratic_zEquation, z);
% solzFun(lambda) = solz;
% negImagIndicies = find(imag(double(solzFun(0.3))) < 0);
% z1 = solz(negImagIndicies(1));
% z2 = solz(negImagIndicies(2));

characteristic_eq = (L2121*L2222*x^2 + (L2121^2 + L2222*L1111 - (L2112 + L2211)^2)*x + L1111*L2121) == 0;
phi_squared = solve(characteristic_eq, x);

z1(lambda) = -1i*simplify(sqrt(-phi_squared(1)));
z2(lambda) = -1i*simplify(sqrt(-phi_squared(2)));

% construct S matrix
S11 = z1*L1212 - L1221*(L1111 + z1^2*L1212)/(z1*(L1122 + L1221));
S12 = z2*L1212 - L1221*(L1111 + z2^2*L1212)/(z2*(L1122 + L1221));
S21 = L2211 - L2222*(L1111 + z1^2*L1212)/(L1122 + L1221);
S22 = L2211 - L2222*(L1111 + z2^2*L1212)/(L1122 + L1221);

solution_eq(lambda) = S11*S22 - S12*S21;
lambdas = 0.01:0.001:0.5;
detS_vect = zeros(size(lambdas));
for k = 1:length(lambdas)
    detS_vect(k) = double(solution_eq(lambdas(k)));
end

plot(lambdas, zeros(size(lambdas)), 'k');
hold on;
plot(lambdas, imag(detS_vect));
ylabel('det(S)    ', 'rot', 1,  'fontsize', 13)
xlabel('\lambda', 'fontsize', 13)
grid on
%%
% % solve for the lambdas that give a nontrivial solution
% solS = double(solve((det(S) == 0), lambda));
% 
% % critical lambda is the lowest positive root:
% 
% solS = real(solS(find(abs(imag(solS)) < 1e-3)));
% lambda_c = min(solS(find(solS > 0)));
% 
% % calculate wave number:
%

wavenumber = 0;
end

