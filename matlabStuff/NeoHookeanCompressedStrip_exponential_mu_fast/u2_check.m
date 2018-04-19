function [] = u2_check(nu, k, w_c, lambda_c, A, alphas, C, x_eval)
mu0 = 1.0;
L = 1;

% give lambda1, the stretch in the x1 direction
lambda1 = 1 - lambda_c;

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


% Now need to make the components of the 6th order M_ijklmn
M = get_d3Wc_dFdFdF(mu0, nu, lambda_c);

% and the v1 and v2 and their derivs
B = zeros(4, 1);
for i = 1:4
    B(i) = (L1111*w_c^2 - k*L1212*alphas(i) - L1212*alphas(i)^2)/(w_c*alphas(i)*(L2112 + L1122) + w_c*k*L2112);
   % B(i) = (w*alphas(i)*(L1221 + L1122) + kappa*L1122*w)/(kappa*L2222*alphas(i) + L2222*alphas(i)^2 - L2121*w^2);
end
v1 = @(x) real(sum(A.*(exp(alphas*x))));
v1_1 = @(x) real(sum(alphas.*(A.*(exp(alphas*x)))));
v1_11 = @(x) real(sum( (alphas.^2).*(A.*(exp(alphas*x)))));
v1_111 = @(x) real(sum( (alphas.^3).*(A.*(exp(alphas*x)))));

v2 = @(x) real(sum((A.*B).*(exp(alphas*x))));
v2_1 = @(x) real(sum( alphas.*((A.*B).*(exp(alphas*x)))));
v2_11 = @(x) real(sum( (alphas.^2).*((A.*B).*(exp(alphas*x)))));
v2_111 = @(x) real(sum( (alphas.^3).*((A.*B).*(exp(alphas*x)))));
% now need to enter all of our functions... (E1, E2, E1', E2', ...)
    function out = E1(x)
        u1 = v1(x);
        u1_1 = v1_1(x);
        u1_11 = v1_11(x);
        u2 = v2(x);
        u2_1 = v2_1(x);
        u2_11 = v2_11(x);
        out =  -M(1,1,1,1,1,1)*(w_c^3)*u1.*u1 + 2*M(1,1,1,1,2,2)*(w_c^2)*u1.*u2_1 - M(1,1,2,2,2,2)*(w_c)*u2_1.*u2_1 + ...
                2*M(1,1,1,2,2,1)*(w_c^2)*u1_1.*u2 + M(1,1,2,1,2,1)*(w_c^3)*u2.*u2 + ...
                k*(M(1,2,1,1,1,2)*w_c*u1.*u1_1 + M(1,2,1,1,2,1)*(w_c^2)*u1.*u2 - M(1,2,2,2,1,2)*u1_1.*u2_1 - M(1,2,2,2,2,1)*w_c*u2.*u2_1) + ... 
                M(1,2,1,1,1,2)*w_c*(u1_1.*u1_1 + u1.*u1_11) + M(1,2,1,1,2,1)*(w_c^2)*(u1_1.*u2 + u1.*u2_1) + ...
                -M(1,2,2,2,1,2)*(u1_11.*u2_1 + u1_1.*u2_11) - M(1,2,2,2,2,1)*w_c*(u2_1.*u2_1 + u2.*u2_11);
    end

    function out = E1_1(x)
        u1 = v1(x);
        u1_1 = v1_1(x);
        u1_11 = v1_11(x);
        u1_111 = v1_111(x);
        u2 = v2(x);
        u2_1 = v2_1(x);
        u2_11 = v2_11(x);
        u2_111 = v2_111(x);

        out = -2*M(1,1,1,1,1,1)*(w_c^3)*u1.*u1_1 + 2*M(1,1,1,1,2,2)*(w_c^2)*(u1_1.*u2_1 + u1.*u2_11) - 2*M(1,1,2,2,2,2)*w_c*u2_1.*u2_11 + ...
            + 2*M(1,1,1,2,1,2)*w_c*u1_1.*u1_11 + 2*M(1,1,1,2,2,1)*(w_c^2)*(u1_11.*u2 + u1_1.*u1_1) + 2*M(1,1,2,1,2,1)*(w_c^3)*u2.*u2_1 + ... 
            k*(M(1,2,1,1,1,2)*w_c*(u1_1.*u1_1 + u1.*u1_11) + M(1,2,1,1,2,1)*(w_c^2)*(u1_1.*u2 + u1.*u2_1) + ...
                -M(1,2,2,2,1,2)*(u1_11.*u2_1 + u1_1.*u2_11) - M(1,2,2,2,2,1)*w_c*(u2_1.*u2_1 + u2.*u2_11)) + ...
            M(1,2,1,1,1,2)*w_c*(u1_1.*u1_11 + u1.*u1_111) + M(1,2,1,1,2,1)*(w_c^2)*(u1_11.*u2 + 2*u1_1.*u2_1 + u1.*u2_11) + ...
            -M(1,2,2,2,1,2)*(u1_111.*u2_1 + 2*u1_11.*u2_11 + u1_1.*u2_111) - M(1,2,2,2,2,1)*w_c*(3*u2_1.*u2_11 + u2.*u2_111);
        
    end

    function out = E2(x)
        u1 = v1(x);
        u1_1 = v1_1(x);
        u1_11 = v1_11(x);
        u2 = v2(x);
        u2_1 = v2_1(x);
        u2_11 = v2_11(x);

        out = 2*(M(2,1,1,1,1,2)*(w_c^2)*u1.*u1_1 + M(2, 1,1,1,2,1)*(w_c^3)*u1.*u2 - M(2,1,2,2,2,2,1)*(w_c^2)*u2.*u2_1 - M(2,1,2,2,1,2)*w_c*u1_1.*u2_1) + ...
            k*(0.5*M(2,2,2,2,2,2)*u2_1.*u2_1 - M(2,2,1,1,2,2)*w_c*u1.*u2_1 + 0.5*M(2,2,1,1,1,1)*(w_c^2)*u1.*u1 - 0.5*M(2,1,2,2,1,2)*u1_1.*u1_1 ... 
                -0.5*M(2,2,2,1,2,1)*(w_c^2)*u2.*u2 - M(2,2,1,2,2,1)*(w_c)*u1_1.*u2) + ...
            M(2,2,2,2,2,2)*u2_1.*u2_11 - M(2,2,1,1,2,2)*w_c*(u1_1.*u2_1 + u1.*u1_11) + M(2,2,1,1,1,1)*(w_c^2)*u1.*u1_1 + ...
            -M(2,2,1,2,1,2)*u1_1.*u1_11 - M(2,2,2,1,2,1)*(w_c^2)*u2.*u2_1 - M(2,2,1,2,2,1)*w_c*(u1_11.*u2 + u1_1.*u2_1);
    end

    function out = E2_1(x)
        u1 = v1(x);
        u1_1 = v1_1(x);
        u1_11 = v1_11(x);
        u1_111 = v1_111(x);
        u2 = v2(x);
        u2_1 = v2_1(x);
        u2_11 = v2_11(x);
        u2_111 = v2_111(x);

        out = 2*(M(2,1,1,1,1,2)*(w_c^2)*(u1_1.*u1_1 + u1.*u1_11) + M(2,1,1,1,2,1)*(w_c^3)*(u1_1.*u2 + u1.*u2_1) - ...
                 M(2,1,2,2,2,1)*(w_c^2)*(u2_1.*u2_1 + u2.*u2_11) - M(2,1,2,2,1,2)*(w_c)*(u1_11.*u2_1 + u1_1.*u2_11)) + ...
              k*(M(2,2,2,2,2,2)*u2_1.*u2_11 - M(2,2,1,1,2,2)*w_c*(u1_1.*u2_1 + u1.*u1_11) + M(2,2,1,1,1,1)*(w_c^2)*u1.*u1_1 + ...
                -M(2,2,1,2,1,2)*u1_1.*u1_11 - M(2,2,2,1,2,1)*(w_c^2)*u2.*u2_1 - M(2,2,1,2,2,1)*w_c*(u1_11.*u2 + u1_1.*u2_1)) + ...
              M(2,2,2,2,2,2)*(u2_11.*u2_11 + u2_1.*u2_111) - M(2,2,1,1,2,2)*(w_c)*(2*u1_1.*u2_11 + u1_11.*u2_1 + u1.*u2_111) + ...
              M(2,2,1,1,1,1)*(w_c^2)*(u1_1.*u1_1 + u1_1.*u1_11) - M(2,2,1,2,1,2)*(u1_11.*u1_11 + u1_1.*u1_111) - ...
              M(2,2,2,1,2,1)*(w_c^2)*(u2_1.*u2_1 + u2.*u2_11) - M(2,2,1,2,2,1)*w_c*(u1_111.*u2 + 2*u1_11.*u2_1 + u1_1.*u2_11);
        
    end

    function out = E2_tilde(x)
        u1 = v1(x);
        u1_1 = v1_1(x);
        u2 = v2(x);
        u2_1 = v2_1(x);
        
        out = 0.5*(M(2,2,2,2,2,2)*u2_1.*u2_1 - 2*M(2,2,1,1,2,2)*w_c*u1.*u2_1 + M(2,2,1,1,1,1)*(w_c^2)*u1.*u1 +...
                   M(2,2,1,2,1,2)*u1_1.*u1_1 + M(2,2,2,1,2,1)*(w_c^2)*u2.*u2 + M(2,2,1,2,2,1)*w_c*u1_1.*u2);
    end

    function out = E2_tilde_1(x)
        u1 = v1(x);
        u1_1 = v1_1(x);
        u1_11 = v1_11(x);
        u2 = v2(x);
        u2_1 = v2_1(x);
        u2_11 = v2_11(x);
        
        out = M(2,2,2,2,2,2)*u2_1.*u2_11 - M(2,2,1,1,2,2)*w_c*(u1_1.*u2_1 + u1.*u2_11) + M(2,2,1,1,1,1)*(w_c^2)*u1.*u1_1 +...
                   M(2,2,1,2,1,2)*u1_1.*u1_11 + M(2,2,2,1,2,1)*(w_c^2)*u2.*u2_1 + M(2,2,1,2,2,1)*w_c*(u1_11.*u2 + u1_1.*u2_1) ;
    end


% now need to enter that A matrix
A_mat = [         0                             1                           0                       0                   ;
    4*(w_c^2)*L1111/L1212                  -k                (2*w_c)*k*L2112/L1212   (2*w_c)*(L2112 + L1122)/L1212  ;
              0                             0                           0                       1                   ;  
    -(2*w_c)*k*L1122/L2222  -(2*w_c)*(L1122 + L2112)/L2222   4*(w_c^2)*L2121/L2222              -k                  ];


[phi, vals] = eigs(A_mat);
phi_inv = inv(phi);
r = diag(vals);

% and define the g's and h's
g = @(x) [zeros(1, length(x)); -E1(x)/L1212; zeros(1, length(x)); -E2(x)/L2222];
h = @(s, x) phi_inv(s, :)*g(x);
g_prime = @(x) [zeros(1, length(x)); -E1_1(x)/L1212; zeros(1, length(x)); -E2_1(x)/L2222];
h_prime = @(s, x) phi_inv(s, :)*g_prime(x);

% now need to setup the w1p and w2p functions
    function out = w1_ph(x)
        out = 0;
        for p = 1:4
            out = out + phi(1, p)*exp(r(p)*x)*integral( @(y) exp(-r(p)*y).*h(p, y), 0, x) + phi(1, p)*C(p)*exp(r(p)*x);
        end
        out = real(out);
    end
    
    function out = w1_ph_1(x)
        out = 0;
        for p = 1:4
            out = out + phi(1, p)*(r(p)*exp(r(p)*x)*integral( @(y) exp(-r(p)*y).*h(p, y), 0, x) + h(p, x)) + phi(1, p)*r(p)*C(p)*exp(r(p)*x);
        end
        out = real(out);
    end

    function out = w1_ph_11(x)
        out = 0;
        for p = 1:4
            out = out + phi(1, p)*(r(p)*r(p)*exp(r(p)*x)*integral( @(y) exp(-r(p)*y).*h(p, y), 0, x) + r(p)*h(p, x) + h_prime(p, x)) + phi(1, p)*r(p)*r(p)*C(p)*exp(r(p)*x);
        end
        out = real(out);
    end

    function out = w2_ph(x)
        out = 0;
        for p = 1:4
            out = out + phi(3, p)*exp(r(p)*x)*integral( @(y) exp(-r(p)*y).*h(p, y), 0, x)+ phi(3, p)*C(p)*exp(r(p)*x) ;
        end
        out = real(out);
    end

    function out = w2_ph_1(x)
        out = 0;
        for p = 1:4
            out = out + phi(3, p)*(r(p)*exp(r(p)*x)*integral( @(y) exp(-r(p)*y).*h(p, y), 0, x) + h(p, x)) + phi(3, p)*r(p)*C(p)*exp(r(p)*x) ;
        end
        out = real(out);
    end

    function out = w2_ph_11(x)
        out = 0;
        for p = 1:4
            out = out + phi(3, p)*(r(p)*r(p)*exp(r(p)*x)*integral( @(y) exp(-r(p)*y).*h(p, y), 0, x) + r(p)*h(p, x) + h_prime(p, x)) + phi(3, p)*r(p)*r(p)*C(p)*exp(r(p)*x) ;
        end
        out = real(out);
    end

%     function out = w2_tilde(x)
%         out = zeros(size(x));
%         for s = 1:length(out)
%             out(s) = -(1/L2222)*integral(@(y) E2_tilde(y), 0, x(s));
%     
%         end
%     end

    function out = w2_tilde_1(x)
        out = -(1/L2222)*E2_tilde(x);
    end

    function out = w2_tilde_11(x)
        out = -(1/L2222)*E2_tilde_1(x);
    end

w1 = zeros(size(x_eval));
w2 = zeros(size(x_eval));
w1_1 = zeros(size(x_eval));
w1_11 = zeros(size(x_eval));
w2_1 = zeros(size(x_eval));
w2_11 = zeros(size(x_eval));

for i = 1:length(x_eval)
   w1(i) = w1_ph(x_eval(i));
   w1_1(i) = w1_ph_1(x_eval(i));
   w1_11(i) = w1_ph_11(x_eval(i));
   w2(i) = w2_ph(x_eval(i));
   w2_1(i) = w2_ph_1(x_eval(i));
   w2_11(i) = w2_ph_11(x_eval(i));
end

LHS1 = -2*w_c*k*L2112*w2 + k*L1212*w1_1 - L1111*4*(w_c^2)*w1 + L1212*w1_11 - (L2112 + L1122)*(2*w_c)*w2_1;
LHS2 = (2*w_c)*k*L1122*w1 + k*L2222*w2_1 + L2222*w2_11 - 4*(w_c^2)*L2121*w2 + (L1122 + L1221)*(2*w_c)*w1_1;
RHS1 = -E1(x_eval);
RHS2 = -E2(x_eval);

figure(1);
hold on;
plot(x_eval, LHS1, 'r', 'Linewidth', 1.2);
grid on
plot(x_eval, RHS1, 'k', 'Linewidth', 1.2);
legend('LHS', 'RHS');
title('Plot of solution for first system of ODEs');
xlabel('x2');


figure(2);
hold on;
plot(x_eval, LHS2, 'r', 'Linewidth', 1.2);
grid on
plot(x_eval, RHS2, 'k', 'Linewidth', 1.2);
legend('LHS', 'RHS');
title('Plot of solution for second system of ODEs');
xlabel('x2');

%% now need to test to see if the complete solution works...

v_ex_1_1 = @(x1, x2) 2*w_c*w1_ph(x2).*cos(2*w_c*x1);
v_ex_1_2 = @(x1, x2) w1_ph_1(x2).*sin(2*w_c*x1);
v_ex_1_12 = @(x1, x2) 2*w_c*w1_ph_1(x2).*cos(2*w_c*x1);
v_ex_1_11 = @(x1, x2) -4*(w_c^2)*w1_ph(x2).*sin(2*w_c*x1);
v_ex_1_22 = @(x1, x2) w1_ph_11(x2).*sin(2*w_c*x1);
v_ex_2_1 = @(x1, x2) -2*w_c*w2_ph(x2).*sin(2*w_c*x1);
v_ex_2_2 = @(x1, x2) w2_ph_1(x2).*cos(2*w_c*x1) + w2_tilde_1(x2);
v_ex_2_12 = @(x1, x2) -2*w_c*w2_ph_1(x2).*sin(2*w_c*x1);
v_ex_2_11 = @(x1, x2) -4*(w_c^2)*w2_ph(x2).*cos(2*w_c*x1);
v_ex_2_22 = @(x1, x2) w2_ph_11(x2).*cos(2*w_c*x1) + w2_tilde_11(x2);

% u1 solution
u1_1_1 = @(x1, x2) -w_c*v1(x2).*cos(w_c*x1);
u1_1_2 = @(x1, x2) -v1_1(x2).*sin(w_c*x1);
u1_1_12 = @(x1, x2) -w_c*v1_1(x2).*cos(w_c*x1);
u1_1_11 = @(x1, x2) (w_c^2)*v1(x2).*sin(w_c*x1);
u1_1_22 = @(x1, x2) -v1_11(x2).*sin(w_c*x1);
u1_2_1 = @(x1, x2) -w_c*v2(x2).*sin(w_c*x1);
u1_2_2 = @(x1, x2) v2_1(x2).*cos(w_c*x1);
u1_2_12 = @(x1, x2) -w_c*v2_1(x2).*sin(w_c*x1);
u1_2_11 = @(x1, x2) -(w_c^2)*v2(x2).*cos(w_c*x1);
u1_2_22 = @(x1, x2) v2_11(x2).*cos(w_c*x1);


[x1_mesh, x2_mesh] = meshgrid([0:0.1:1], [0:0.01:1]);
LHS1_PDE = zeros(size(x1_mesh));
LHS2_PDE = zeros(size(x1_mesh));
RHS1_PDE = zeros(size(x1_mesh));
RHS2_PDE = zeros(size(x1_mesh));

lhs_pde_1 = @(x1, x2) k*L2112*v_ex_2_1(x1, x2) + k*L1212*v_ex_1_2(x1, x2) + ...
                      L1111*v_ex_1_11(x1, x2) + L1212*v_ex_1_22(x1, x2) + ...
                      (L2112 + L1122)*v_ex_2_12(x1, x2);
               
lhs_pde_2 = @(x1, x2) k*L1122*v_ex_1_1(x1, x2) + k*L2222*v_ex_2_2(x1, x2) + ...
                      L2222*v_ex_2_22(x1, x2) + L1212*v_ex_2_11(x1, x2) + ...
                      (L2112 + L1122)*v_ex_1_12(x1, x2);
                  
    function out = rhs_pde_1(x1, x2)
        a1_1 = u1_1_1(x1, x2);
        a1_2 = u1_1_2(x1, x2);
        a1_12 = u1_1_12(x1, x2);
        a1_11 = u1_1_11(x1, x2);
        a1_22 = u1_1_22(x1, x2);
        a2_1 = u1_2_1(x1, x2);
        a2_2 = u1_2_2(x1, x2);
        a2_12 = u1_2_12(x1, x2);
        a2_11 = u1_2_11(x1, x2);
        a2_22 = u1_2_22(x1, x2);

        out = -2*(M(1,1,1,1,1,1)*a1_1*a1_11 + M(1,1,1,1,2,2)*(a1_11*a2_2 + a1_1*a2_12) + ...
                  M(1,1,2,2,2,2)*a2_2*a2_12 + M(1,1,1,2,1,2)*a1_2*a1_12 + M(1,1,1,2,2,1)*(a1_12*a2_1 + a1_2*a2_11) + ...
                  M(1,1,2,1,2,1)*a2_1*a2_11 + ...
                  k*(M(1,2,1,1,1,2)*a1_1*a1_2 + M(1,2,1,1,2,1)*a1_1*a2_1 + M(1,2,2,2,1,2)*a2_2*a1_2 + M(1,2,2,2,2,1)*a2_2*a2_1) + ...
                  (M(1,2,1,1,1,2)*(a1_12*a1_2 + a1_1*a1_22) + M(1,2,1,1,2,1)*(a1_12*a2_1 + a1_1*a2_12) + ...
                   M(1,2,2,2,1,2)*(a2_22*a1_2 + a2_2*a1_22) + M(1,2,2,2,2,1)*(a2_22*a2_1 + a2_2*a2_12)));
    end

    function out = rhs_pde_2(x1, x2)
        a1_1 = u1_1_1(x1, x2);
        a1_2 = u1_1_2(x1, x2);
        a1_12 = u1_1_12(x1, x2);
        a1_11 = u1_1_11(x1, x2);
        a1_22 = u1_1_22(x1, x2);
        a2_1 = u1_2_1(x1, x2);
        a2_2 = u1_2_2(x1, x2);
        a2_12 = u1_2_12(x1, x2);
        a2_11 = u1_2_11(x1, x2);
        a2_22 = u1_2_22(x1, x2);
        
        out = -(2*(M(2,1,1,1,1,2)*(a1_11*a1_2 + a1_1*a1_12) + M(2,1,1,1,2,1)*(a1_11*a2_1 + a1_1*a2_11) + ...
                  M(2,1,2,2,2,1)*(a2_12*a2_1 + a2_2*a2_11) + M(2,1,2,2,1,2)*(a2_12*a1_2 + a2_2*a1_12)) + ...
                k*(M(2,2,2,2,2,2)*a2_2*a2_2 + 2*M(2,2,1,1,2,2)*a1_1*a2_2 + M(2,2,1,1,1,1)*a1_1*a1_1 + ...
                   M(2,2,1,2,1,2)*a1_2*a1_2 + M(2,2,2,1,2,1)*a2_1*a2_1 + 2*M(2,2,1,2,2,1)*a1_2*a2_1) + ...
                2*(M(2,2,2,2,2,2)*a2_2*a2_22 + M(2,2,1,1,2,2)*(a1_12*a2_2 + a1_1*a2_22) + ...
                   M(2,2,1,1,1,1)*a1_1*a1_12 + M(2,2,1,2,1,2)*a1_2*a1_22 + + ...
                   M(2,2,2,1,2,1)*a2_1*a2_12 + M(2,2,1,2,2,1)*(a1_22*a2_1 + a1_2*a2_12)));
    end

for i = 1:length(x1_mesh)
    for j = 1:length(x2_mesh)
        LHS1_PDE(i, j) = lhs_pde_1(x1_mesh(i, j), x2_mesh(i, j));
        LHS2_PDE(i, j) = lhs_pde_2(x1_mesh(i, j), x2_mesh(i, j));
        RHS1_PDE(i, j) = rhs_pde_1(x1_mesh(i, j), x2_mesh(i, j));
        RHS2_PDE(i, j) = rhs_pde_2(x1_mesh(i, j), x2_mesh(i, j));
    end
end


figure(3) 
surf(x1_mesh, x2_mesh, LHS1_PDE);
hold on;
surf(x1_mesh, x2_mesh, RHS1_PDE);
grid on;
xlabel('x1');
ylabel('x2');
legend('lhs', 'rhs');
title('lhs and rhs for pde 1');

figure(4) 
surf(x1_mesh, x2_mesh, LHS2_PDE);
hold on;
surf(x1_mesh, x2_mesh, RHS2_PDE);
grid on;
xlabel('x1');
ylabel('x2');
legend('lhs', 'rhs');
title('lhs and rhs for pde 2');

end