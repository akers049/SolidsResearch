function [C] = u2_exponential_mu(nu, k, w_c, lambda_c, A, alphas)
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
v2 = @(x) real(sum((A.*B).*(exp(alphas*x))));
v2_1 = @(x) real(sum( alphas.*((A.*B).*(exp(alphas*x)))));
v2_11 = @(x) real(sum( (alphas.^2).*((A.*B).*(exp(alphas*x)))));

% now need to enter all of our functions... (E1, E2, F1, F2)
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

    function out = F1(x)
        u1 = v1(x);
        u1_1 = v1_1(x);
        u2 = v2(x);
        u2_1 = v2_1(x);
        
        out = M(1,2,1,1,1,2)*w_c*u1*u1_1 + M(1,2,1,1,2,1)*(w_c^2)*u1*u2 - M(1,2,2,2,1,2)*u1_1*u2_1 - M(1,2,2,2,2,1)*w_c*u2*u2_1;
    end

    function out = F2(x)
        u1 = v1(x);
        u1_1 = v1_1(x);
        u2 = v2(x);
        u2_1 = v2_1(x);
        
        out = 0.5*(M(2,2,2,2,2,2)*u2_1*u2_1 - 2*M(2,2,1,1,2,2)*w_c*u1*u2_1 + M(2,2,1,1,1,1)*(w_c^2)*u1*u1 + ...
                   -M(2,2,1,2,1,2)*u1_1*u1_1 - M(2,2,2,1,2,1)*(w_c^2)*u2*u2 - M(2,2,1,2,2,1)*(w_c)*u1_1*u2);
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

% now need to setup the w1p and w2p functions
    function out = w1p(x)
        out = 0;
        for p = 1:4
            out = out + phi(1, p)*exp(r(p)*x)*integral( @(y) exp(-r(p)*y).*h(p, y), 0, x);
        end
        out = real(out);
    end

    function out = w1p_1(x)
        out = 0;
        for p = 1:4
            out = out + phi(1, p)*(r(p)*exp(r(p)*x)*integral( @(y) exp(-r(p)*y).*h(p, y), 0, x) + h(p, x));
        end
        out = real(out);
    end

    function out = w2p(x)
        out = 0;
        for p = 1:4
            out = out + phi(3, p)*exp(r(p)*x)*integral( @(y) exp(-r(p)*y).*h(p, y), 0, x);
        end
        out = real(out);
    end

    function out = w2p_1(x)
        out = 0;
        for p = 1:4
            out = out + phi(3, p)*(r(p)*exp(r(p)*x)*integral( @(y) exp(-r(p)*y).*h(p, y), 0, x) + h(p, x));
        end
        out = real(out);
    end


% and now we are finally ready to assemble the boundary condition matrix
% and the rhs
boundary_mat = zeros(4, 4);
for i = 1:4
   boundary_mat(1, i) = phi(3, i);
   boundary_mat(2, i) = phi(1, i)*L1212*r(i) - phi(3, i)*L1221*2*w_c;
   boundary_mat(3, i) = phi(1, i)*L1212*r(i)*exp(r(i)*L) - phi(3, i)*L1221*(2*w_c)*exp(r(i)*L);
   boundary_mat(4, i) = phi(1, i)*L1122*(2*w_c)*exp(r(i)*L) + phi(3, i)*L2222*r(i)*exp(r(i)*L);
end
rhs = [ 0    ;
        -F1(0);
        -F1(L) - L1212*w1p_1(L) + L1221*(2*w_c)*w2p(L) ;
        -F2(L) - L1122*(2*w_c)*w1p(L) - L2222*w2p_1(L) ];
    
    
C = boundary_mat\rhs;
end