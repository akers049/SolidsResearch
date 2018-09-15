function [PD] = do_asymptotics_exponential_and_print_file(nu, kappa, starting_wavelength, index)


PD = struct;

kappa_eval = kappa;
stepSize = 0.05;
wavelength_vector = starting_wavelength:stepSize:60;
solutionVector = zeros(size(wavelength_vector));

% get initial guess
xx = 0.1:0.001:0.89;
yy = zeros(size(xx));
for i = 1:length(xx)
    yy(i) = evaluate_determinant_M(wavelength_vector(1), nu, kappa_eval, xx(i));
end
% figure(1);
% plot(xx, yy, 'LineWidth',1.3);
% axis([0.1 0.99 -min(abs(yy)) min(abs(yy))]);
% grid on;
% [guess_x, ~] = ginput(1);
% close all;
yy_sign = sign(yy);
crossingIndex = find(diff(sign(yy_sign), 1));
guess_x = xx(crossingIndex(end));


MAXCOUNT = 30;
initialGuess = guess_x;

ds = stepSize/100;
for i = 1:length(wavelength_vector)
    funHandle = @(x) evaluate_determinant_M(wavelength_vector(i), nu, kappa_eval, x);

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
%       fprintf('%f NEXT\n', i);
      solutionVector(i) = initialGuess;
    end
%     guess_low = initialGuess - stepSize/5;
%     guess_high = initialGuess + stepSize/5; 
%     if(sign(funHandle(guess_low)) == sign(funHandle(guess_high)))
%         solutionVector(i) = initialGuess;
%     else
%         solutionVector(i) = fzero(funHandle, [guess_low, guess_high]);
%     end
    
    if ( i >1)
        initialGuess = 2*solutionVector(i) - solutionVector(i-1);
    else
        initialGuess = solutionVector(i);
    end
    
    if(solutionVector(i) > 0.8)
        break;
    end
    
end
solutionVector = solutionVector(1:i);
wavelength_vector = wavelength_vector(1:i);

plot(wavelength_vector, solutionVector, 'LineWidth', 1.3)
grid on
ylabel('$\lambda_{c}$      ', 'rot', 1, 'Interpreter', 'Latex', 'Fontsize', 16)
xlabel('Wavelength', 'fontsize', 16)
title(['$\kappa$ = ', num2str(kappa)], 'Fontsize', 16, 'Interpreter', 'Latex')
axis([0 35 0.2 1])

criticalLoad = min(solutionVector);
criticalWavelength =  wavelength_vector(find(min(solutionVector) == solutionVector));
criticalFreq = 2*pi/criticalWavelength;

if(criticalLoad < 0 || criticalLoad > 1)
   PD.good = 0;
   return; 
end

%% now get the A and alphas
lambda_eval = criticalLoad;
wavelength_eval = criticalWavelength;
kappa_eval = kappa;

w_eval = 2*pi/wavelength_eval;

mu0 = 1;
L = 1;

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
b = (2*kappa_eval*L2222*L1212);
c = (w_eval^2*(L1122^2 + 2*L1122*L2112 + L2112^2 - L1212^2  - L1111*L2222) + L2222*L1212*(kappa_eval^2));
d = (kappa_eval*w_eval^2*(L1122^2 + 2*L2112*L1122  + L2112^2 - L1212^2 - L1111*L2222));
e = (w_eval^2*(L1122*L2112*kappa_eval^2 + L1111*L1212*w_eval^2));

% solve for alphas
alphas_roots = roots([a, b, c, d, e]);
alphas = sort(alphas_roots);
% if (any(abs(imag(alphas)) > 0))
%    detVal = 1;
%    return;
% end
alphas
% now can solve for B's:
B = zeros(4, 1);
for i = 1:4
    B(i) = (L1111*w_eval^2 - kappa_eval*L1212*alphas(i) - L1212*alphas(i)^2)/(w_eval*alphas(i)*(L2112 + L1122) + w_eval*kappa_eval*L2112);
   % B(i) = (w*alphas(i)*(L1221 + L1122) + kappa*L1122*w)/(kappa*L2222*alphas(i) + L2222*alphas(i)^2 - L2121*w^2);
end

% now construct the system matrix:

systemMat = (zeros(4,4));
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
    systemMat(3, i) = L2222*alphas(i)*B(i)*exp(alphas(i)*L) - w_eval*L1122*exp(alphas(i)*L);
end

% fourth row
for i = 1:4
    systemMat(4, i) = w_eval*L2112*B(i)*exp(alphas(i)*L) + L1212*alphas(i)*exp(alphas(i)*L);
end

systemMat_reffed = rref(systemMat);
A = [systemMat_reffed(1:3, 4); -1];
if (abs(imag(A(1))) > 1e-6)
    val = imag(A(2)) - imag(A(1)) - sqrt(-1)*(real(A(1,1)) - real(A(2)));
    A = A*val;
end
A = A/sum(B.*A(1:4).*exp(alphas));

norm_2_v1_v2 = @(x) (sum(A(1:4).*exp(alphas.*x))).^2 + (sum(B.*A(1:4).*exp(alphas.*x))).^2;
u1_dot_u1 = (1/criticalWavelength)*pi/(criticalFreq)*integral(norm_2_v1_v2, 0, 1);
% u1_dot_u1 = pi/(criticalFreq)*integral(norm_2_v1_v2, 0, 1);

A = A*sqrt(1/u1_dot_u1);

PD.k = kappa;
PD.lambda_c = criticalLoad;
PD.w_c = criticalFreq;
PD.A = A;
PD.B = B;
PD.alphas = alphas;
PD.nu = nu;
PD.A = A;
% u1u1 = @(x1, x2) testFunc(PD, x1, x2);
% u1_dot_u1 = integral2(u1u1, -pi/PD.w_c, pi/PD.w_c, 0, 1);
% 
% A = A*sqrt(1/u1_dot_u1);
% PD.A = A;


PD = u2_exponential_mu(PD);

%% print the file

o = fopen(['/home/andrew/dealii/SolidsResearch/CompressedStrip_Project/asymptoticCalculation/exp_double_root_finding/', 'asymptotic_exponential_', num2str(index), '.in'],'w');
fprintf(o, '#INPUT FILE FOR ASYMPTOTIC EXPONENTIAL mu\n');
fprintf(o, 'exponential_mu\n');
fprintf(o, '1\n');
fprintf(o, '15 15\n');
fprintf(o, '%f %f\n', PD.nu, PD.k);
fprintf(o, '%f\n', PD.lambda_c);
fprintf(o, '%f\n', PD.w_c);
for i = 1:4
   fprintf(o, '%.16f %.16f ', real(PD.alphas(i)), imag(PD.alphas(i)));  
end
fprintf(o, '\n');

for i = 1:4
   fprintf(o, '%.16f %.16f ', real(PD.A(i)), imag(PD.A(i)));  
end
fprintf(o, '\n');

for i = 1:4
   fprintf(o, '%.16f %.16f ', real(PD.B(i)), imag(PD.B(i)));  
end
fprintf(o, '\n');

for i = 1:4
   fprintf(o, '%.16f %.16f ', real(PD.phi1(i)), imag(PD.phi1(i)));  
end
fprintf(o, '\n');

for i = 1:4
   fprintf(o, '%.16f %.16f ', real(PD.phi3(i)), imag(PD.phi3(i)));  
end
fprintf(o, '\n');

for i = 1:4
   fprintf(o, '%.16f %.16f ', real(PD.phi_inv_T_2(i)), imag(PD.phi_inv_T_2(i)));  
end
fprintf(o, '\n');

for i = 1:4
   fprintf(o, '%.16f %.16f ', real(PD.phi_inv_T_4(i)), imag(PD.phi_inv_T_4(i)));  
end
fprintf(o, '\n');

for i = 1:4
   fprintf(o, '%.16f %.16f ', real(PD.r(i)), imag(PD.r(i)));  
end
fprintf(o, '\n');

for i = 1:4
   fprintf(o, '%.16f %.16f ', real(PD.C(i)), imag(PD.C(i)));  
end
fprintf(o, '\n');

for i = 1:4
   fprintf(o, '%f %f ', real(PD.C(i)), imag(PD.C(i)));  
end
fprintf(o, '\n');

fclose(o);

PD.good = 1;
end