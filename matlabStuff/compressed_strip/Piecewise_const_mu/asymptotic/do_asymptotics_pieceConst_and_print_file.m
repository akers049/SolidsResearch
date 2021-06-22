function [PD] = do_asymptotics_pieceConst_and_print_file(nu, L1, kappa, starting_wavelength, index, just_critload)

PD = struct;
kappa_eval = kappa;
stepSize = 0.0001;
wavelength_vector = starting_wavelength:stepSize:6;
solutionVector = zeros(size(wavelength_vector));

% get initial guess
xx = 0.002:0.01:0.89;
yy = zeros(size(xx));
for i = 1:length(xx)
    yy(i) = piece_const_mu_evaluate_determinant_M(nu, L1, kappa_eval, wavelength_vector(1), xx(i));
end
% figure(1);
% plot(xx, yy, 'LineWidth',1.3);
% axis([0.1 0.99 -min(abs(yy)) min(abs(yy))]);
% grid on;
% [guess_x, ~] = ginput(1);
% close all;
yy_sign = sign(yy);
crossingIndex = find(diff(sign(yy_sign), 1));
guess_x = xx(crossingIndex(1));

MAXCOUNT = 100;
initialGuess = guess_x;

ds = stepSize;
for i = 1:length(wavelength_vector)
%     funHandle = @(x) piece_const_mu_evaluate_determinant_M(nu, L1, kappa_eval, wavelength_vector(i), x);
% 
%     solutionVector(i) = fzero(funHandle, [initialGuess - 0.1, initialGuess + 0.1]);
%     
%     initialGuess = solutionVector(i);
%     if(solutionVector(i) > 0.8)
%         break;
%     end

    funHandle = @(x) piece_const_mu_evaluate_determinant_M(nu, L1, kappa_eval, wavelength_vector(i), x);
    guess_low = initialGuess - ds;
    guess_high = initialGuess + ds;

    
    count = 0;
    while((sign(funHandle(guess_low)) == sign(funHandle(guess_high))) && count < MAXCOUNT)
        guess_low = guess_low - ds;
        guess_high = guess_high + ds;
        count = count + 1;
    end
        
    
    if(count < MAXCOUNT)
      solutionVector(i) = fzero(funHandle, [guess_low, guess_high]);
    else
       fprintf('%f NEXT\n', i);
      solutionVector(i) = initialGuess;
    end
    
    if ( i >1)
        initialGuess = 2*solutionVector(i) - solutionVector(i-1);
    else
        initialGuess = solutionVector(i);
    end
    
    if(solutionVector(i) > 0.9)
        break;
    end
    
end
solutionVector = solutionVector(1:i);
wavelength_vector = wavelength_vector(1:i);

% plot(wavelength_vector, solutionVector, 'LineWidth', 1.3)
% grid on
% ylabel('$\lambda_{c}$      ', 'rot', 1, 'Interpreter', 'Latex', 'Fontsize', 16)
% xlabel('Wavelength', 'fontsize', 16)
% title(['$\mu_2/\mu_1$ = ', num2str(kappa)], 'Fontsize', 16, 'Interpreter', 'Latex')
% axis([0 25 0 1])


criticalLoad = min(solutionVector);
criticalWavelength =  wavelength_vector(find(min(solutionVector) == solutionVector));
criticalWavelength = criticalWavelength(1);
criticalFreq = 2*pi/criticalWavelength;

if(criticalLoad < 0 || criticalLoad > 1)
   PD.good = 0;
   return; 
end

PD.wavelength_c = criticalWavelength;
PD.lambda_c = criticalLoad;

if(just_critload == true)
    PD.good = 1;
    return;
end
%% Now get the eigenvectors and stuff:

lambda_eval = criticalLoad;
kappa_eval = kappa;
wavelength_eval = criticalWavelength;


if abs(lambda_eval - 1) < 1e-1
    detVal = 1e9;
    return;
end

w_eval = 2*pi/wavelength_eval;
L2 = 1.0;
mu0 = 1;

% give lambda1, the stretch in the x1 direction
lambda1 = 1 - lambda_eval;

% solve for lambda2, the stretch in the x2 direction (use the positive
% solution).
lambda2_sols = roots([(mu0 + lambda1^2*2*nu*mu0/(1 - nu)), -(2*mu0*nu/(1 - nu))*lambda1, -mu0]);
lambda2 = (lambda2_sols(find(double(lambda2_sols >= 0))));

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

% Make the characteristic equation;
a = (L2222*L1212);
c = w_eval^2*(L1122^2 + 2*L1122*L2112 + L2112^2 - L1212^2  - L1111*L2222);
e = w_eval^4*(L1111*L1212);

% solve for alphas
alphas_roots = roots([a, 0, c, 0, e]);
alphas = sort(alphas_roots);

% now can solve for B's:
B = zeros(4, 1);
for i = 1:4
    B(i) = (L1111*w_eval^2 - L1212*alphas(i)^2)/(w_eval*alphas(i)*(L2112 + L1122));
   % B(i) = (w*alphas(i)*(L1221 + L1122) + kappa*L1122*w)/(kappa*L2222*alphas(i) + L2222*alphas(i)^2 - L2121*w^2);
end

% now construct the system matrix:

systemMat = (zeros(8,8));
% first row
for i = 1:4
    systemMat(1, i) = B(i);
end

% second row
for i = 1:4
    systemMat(2, i) = L1212*alphas(i) + w_eval*L2112*B(i);
end

% third row
for i = 1:4
    systemMat(3, i+4) = L2222*alphas(i)*B(i)*exp(alphas(i)*L2) - w_eval*L1122*exp(alphas(i)*L2);
end

% fourth row
for i = 1:4
    systemMat(4, i+4) = w_eval*L2112*B(i)*exp(alphas(i)*L2) + L1212*alphas(i)*exp(alphas(i)*L2);
end

% fifth row
for i = 1:4
   systemMat(5, i) = exp(alphas(i)*L1); 
   systemMat(5, i+4) = -exp(alphas(i)*L1); 
end

% sixth row
for i = 1:4
   systemMat(6, i) = B(i)*exp(alphas(i)*L1); 
   systemMat(6, i+4) = -B(i)*exp(alphas(i)*L1); 
end

% seventh row
for i = 1:4
    systemMat(7, i) = w_eval*L2112*B(i)*exp(alphas(i)*L1) + L1212*alphas(i)*exp(alphas(i)*L1);
    systemMat(7, i+4) = -kappa_eval*(w_eval*L2112*B(i)*exp(alphas(i)*L1) + L1212*alphas(i)*exp(alphas(i)*L1));
end

% eighth row
for i = 1:4
    systemMat(8, i) = L2222*alphas(i)*B(i)*exp(alphas(i)*L1) - w_eval*L1122*exp(alphas(i)*L1);
    systemMat(8, i+4) = -kappa_eval*(L2222*alphas(i)*B(i)*exp(alphas(i)*L1) - w_eval*L1122*exp(alphas(i)*L1));
    
end
det(systemMat);
[eigVects, eigVals] = eig(systemMat);
[~, mindex] = min(abs(diag(eigVals)));
A = eigVects(:, mindex);
% A = null(systemMat);
% A = A(:, 1);
% systemMat_reffed = rref(systemMat);
% A = [systemMat_reffed(1:7, 8); -1];
% A = A/sum(B.*A(5:8).*exp(alphas));

% norm_2_v1_v2_1 = @(x) (sum(A(1:4).*exp(alphas.*x))).^2 + (sum(B.*A(1:4).*exp(alphas.*x))).^2;
% norm_2_v1_v2_2 = @(x) (sum(A(5:8).*exp(alphas.*x))).^2 + (sum(B.*A(5:8).*exp(alphas.*x))).^2;
% 
% u1_dot_u1 = (1/criticalWavelength)*pi/(criticalFreq)*(integral(norm_2_v1_v2_1, 0, L1) + integral(norm_2_v1_v2_2, L1, 1)) ;
% % u1_dot_u1 = pi/(criticalFreq)*integral(norm_2_v1_v2, 0, 1);
% 
% A = A*sqrt(1/u1_dot_u1);

PD.k = kappa;
PD.w_c = criticalFreq;
PD.nu = nu;
PD.L1 = L1;

if ( abs(det(systemMat)) < 1.0e-4)
    PD.A = A;
    PD.B = B;
    PD.alphas = alphas;
    u1_inner = @(x1, x2) testFunc_1(PD, x1, x2);
    u1_dot_u1 = (1/criticalWavelength)*(integral2(u1_inner, -pi/PD.w_c, pi/PD.w_c, 0 , L1) + integral2(u1_inner, -pi/PD.w_c, pi/PD.w_c, L1, 1));
    PD.A = PD.A*sqrt(1/u1_dot_u1);
else
    PD = vpa_zero_and_A(PD, criticalLoad);
    u1_inner = @(x1, x2) testFunc_1(PD, x1, x2);
    u1_dot_u1 = (1/criticalWavelength)*(integral2(u1_inner, -pi/PD.w_c, pi/PD.w_c, 0 , L1) + integral2(u1_inner, -pi/PD.w_c, pi/PD.w_c, L1, 1));
    PD.A = PD.A*sqrt(1/u1_dot_u1);
end
% u1u1 = @(x1, x2) u1_dot_u1(PD, x1, x2);
% u1_dot = integral2(u1u1, -pi/PD.w_c, pi/PD.w_c, 0, 1);
% 
% A = A*sqrt(1/u1_dot);

% PD.A = A;


PD = u2_piece_const_mu(PD);

%% print the file

o = fopen(['/home/andrew/Research/MinnesotaStuff/SolidsResearch/CompressedStrip_Project/asymptoticCalculation/inputFiles/piecewiseConstInputFiles/', 'asymptotic_piecewiseConstant_', num2str(index), '.in'],'w');
fprintf(o, '#INPUT FILE FOR PIECEWISE CONSTANT mu\n');
fprintf(o, 'piecewiseConstant_mu\n');
fprintf(o, '1\n');
fprintf(o, '20 40 1\n');
fprintf(o, '%e %e %e\n', PD.nu, PD.k, PD.L1);
fprintf(o, '%.16e\n', PD.lambda_c);
fprintf(o, '%.16e\n', PD.w_c);
for i = 1:4
   fprintf(o, '%.16e %.16e ', real(PD.alphas(i)), imag(PD.alphas(i)));  
end
fprintf(o, '\n');

for i = 1:8
   fprintf(o, '%.16e %.16e ', real(PD.A(i)), imag(PD.A(i)));  
end
fprintf(o, '\n');

for i = 1:4
   fprintf(o, '%.16e %.16e ', real(PD.B(i)), imag(PD.B(i)));  
end
fprintf(o, '\n');

for i = 1:4
   fprintf(o, '%.16e %.16e ', real(PD.phi1(i)), imag(PD.phi1(i)));  
end
fprintf(o, '\n');

for i = 1:4
   fprintf(o, '%.16e %.16e ', real(PD.phi3(i)), imag(PD.phi3(i)));  
end
fprintf(o, '\n');

for i = 1:4
   fprintf(o, '%.16e %.16e ', real(PD.phi_inv_T_2(i)), imag(PD.phi_inv_T_2(i)));  
end
fprintf(o, '\n');

for i = 1:4
   fprintf(o, '%.16e %.16e ', real(PD.phi_inv_T_4(i)), imag(PD.phi_inv_T_4(i)));  
end
fprintf(o, '\n');

for i = 1:4
   fprintf(o, '%.16e %.16e ', real(PD.r(i)), imag(PD.r(i)));  
end
fprintf(o, '\n');

for i = 1:8
   fprintf(o, '%.16e %.16e ', real(PD.C(i)), imag(PD.C(i)));  
end
fprintf(o, '\n');

fclose(o);

PD.good = 1;

end