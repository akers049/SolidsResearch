function [PD] = u2_piece_const_mu(PD)
nu = PD.nu;
L1 = PD.L1;
k = PD.k;
w_c = PD.w_c;
lambda_c = PD.lambda_c;
A1 = PD.A(1:4);
A2 = PD.A(5:end);
alphas = PD.alphas;

L2 = 1;
mu0 = 1;

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
    B(i) = (L1111*w_c^2 - L1212*alphas(i)^2)/(w_c*alphas(i)*(L2112 + L1122));
   % B(i) = (w*alphas(i)*(L1221 + L1122) + kappa*L1122*w)/(kappa*L2222*alphas(i) + L2222*alphas(i)^2 - L2121*w^2);
end
v1 = @(x, A) real(sum(A.*(exp(alphas*x))));
v1_1 = @(x, A) real(sum(alphas.*(A.*(exp(alphas*x)))));
v1_11 = @(x, A) real(sum( (alphas.^2).*(A.*(exp(alphas*x)))));
v2 = @(x, A) real(sum((A.*B).*(exp(alphas*x))));
v2_1 = @(x, A) real(sum( alphas.*((A.*B).*(exp(alphas*x)))));
v2_11 = @(x, A) real(sum( (alphas.^2).*((A.*B).*(exp(alphas*x)))));


% now need to enter all of our functions... (E1, E2, F1, F2)
    function out = E1(x, dom)
        if(dom == 1)
            u1 = v1(x, A1);
            u1_1 = v1_1(x, A1);
            u1_11 = v1_11(x, A1);
            u2 = v2(x, A1);
            u2_1 = v2_1(x, A1);
            u2_11 = v2_11(x, A1);
        else
            u1 = v1(x, A2);
            u1_1 = v1_1(x, A2);
            u1_11 = v1_11(x, A2);
            u2 = v2(x, A2);
            u2_1 = v2_1(x, A2);
            u2_11 = v2_11(x, A2);
        end
        
        out =  -M(1,1,1,1,1,1)*(w_c^3)*u1.*u1 + 2*M(1,1,1,1,2,2)*(w_c^2)*u1.*u2_1 - M(1,1,2,2,2,2)*(w_c)*u2_1.*u2_1 + ...
                2*M(1,1,1,2,2,1)*(w_c^2)*u1_1.*u2 + M(1,1,2,1,2,1)*(w_c^3)*u2.*u2 + ...
                M(1,2,1,1,1,2)*w_c*(u1_1.*u1_1 + u1.*u1_11) + M(1,2,1,1,2,1)*(w_c^2)*(u1_1.*u2 + u1.*u2_1) + ...
                -M(1,2,2,2,1,2)*(u1_11.*u2_1 + u1_1.*u2_11) - M(1,2,2,2,2,1)*w_c*(u2_1.*u2_1 + u2.*u2_11);
    end

    function out = E2(x, dom)
        if(dom == 1)
            u1 = v1(x, A1);
            u1_1 = v1_1(x, A1);
            u1_11 = v1_11(x, A1);
            u2 = v2(x, A1);
            u2_1 = v2_1(x, A1);
            u2_11 = v2_11(x, A1);
        else
            u1 = v1(x, A2);
            u1_1 = v1_1(x, A2);
            u1_11 = v1_11(x, A2);
            u2 = v2(x, A2);
            u2_1 = v2_1(x, A2);
            u2_11 = v2_11(x, A2);
        end
        
        out = 2*(M(2,1,1,1,1,2)*(w_c^2)*u1.*u1_1 + M(2,1,1,1,2,1)*(w_c^3)*u1.*u2 - M(2,1,2,2,2,1)*(w_c^2)*u2.*u2_1 - M(2,1,2,2,1,2)*w_c*u1_1.*u2_1) + ...
            M(2,2,2,2,2,2)*u2_1.*u2_11 - M(2,2,1,1,2,2)*w_c*(u1_1.*u2_1 + u1.*u2_11) + M(2,2,1,1,1,1)*(w_c^2)*u1.*u1_1 + ...
            -M(2,2,1,2,1,2)*u1_1.*u1_11 - M(2,2,2,1,2,1)*(w_c^2)*u2.*u2_1 - M(2,2,1,2,2,1)*w_c*(u1_11.*u2 + u1_1.*u2_1);
    end

    function out = F1(x, dom)
        if(dom == 1)
            u1 = v1(x, A1);
            u1_1 = v1_1(x, A1);
            u2 = v2(x, A1);
            u2_1 = v2_1(x, A1);
        else
            u1 = v1(x, A2);
            u1_1 = v1_1(x, A2);
            u2 = v2(x, A2);
            u2_1 = v2_1(x, A2);
        end
        
        out = M(1,2,1,1,1,2)*w_c*u1*u1_1 + M(1,2,1,1,2,1)*(w_c^2)*u1*u2 - M(1,2,2,2,1,2)*u1_1*u2_1 - M(1,2,2,2,2,1)*w_c*u2*u2_1;
    end

    function out = F2(x, dom)
        if(dom == 1)
            u1 = v1(x, A1);
            u1_1 = v1_1(x, A1);
            u2 = v2(x, A1);
            u2_1 = v2_1(x, A1);
        else
            u1 = v1(x, A2);
            u1_1 = v1_1(x, A2);
            u2 = v2(x, A2);
            u2_1 = v2_1(x, A2);
        end
        
        out = 0.5*(M(2,2,2,2,2,2)*u2_1*u2_1 - 2*M(2,2,1,1,2,2)*w_c*u1*u2_1 + M(2,2,1,1,1,1)*(w_c^2)*u1*u1 + ...
                   -M(2,2,1,2,1,2)*u1_1*u1_1 - M(2,2,2,1,2,1)*(w_c^2)*u2*u2 - 2*M(2,2,1,2,2,1)*(w_c)*u1_1*u2);
    end

% now need to enter that A matrix
A_mat = [         0                             1                           0                       0                   ;
    4*(w_c^2)*L1111/L1212                  0                          0               (2*w_c)*(L2112 + L1122)/L1212  ;
              0                             0                           0                       1                   ;  
             0             -(2*w_c)*(L1122 + L2112)/L2222   4*(w_c^2)*L2121/L2222              0                  ];


[phi, vals] = eigs(A_mat);
phi_inv = inv(phi);
r = diag(vals);

% and define the g's and h's
g = @(x, dom) [zeros(1, length(x)); -E1(x, dom)/L1212; zeros(1, length(x)); -E2(x, dom)/L2222];
h = @(s, x, dom) phi_inv(s, :)*g(x, dom);

% now need to setup the w1p and w2p functions
    function out = w1p(x, dom)
        out = 0;
        low_bound = L1*(dom -1);
        for p = 1:4
            out = out + phi(1, p)*exp(r(p)*x)*integral( @(y) exp(-r(p)*y).*h(p, y, dom), low_bound, x);
        end
        out = real(out);
    end

    function out = w1p_1(x, dom)
        out = 0;
        low_bound = L1*(dom -1);
        for p = 1:4
            out = out + phi(1, p)*(r(p)*exp(r(p)*x)*integral( @(y) exp(-r(p)*y).*h(p, y, dom), low_bound, x) + h(p, x, dom));
        end
        out = real(out);
    end

    function out = w2p(x, dom)
        out = 0;
        low_bound = L1*(dom -1);
        for p = 1:4
            out = out + phi(3, p)*exp(r(p)*x)*integral( @(y) exp(-r(p)*y).*h(p, y, dom), low_bound, x);
        end
        out = real(out);
    end

    function out = w2p_1(x, dom)
        out = 0;
        low_bound = L1*(dom -1);
        for p = 1:4
            out = out + phi(3, p)*(r(p)*exp(r(p)*x)*integral( @(y) exp(-r(p)*y).*h(p, y, dom), low_bound, x) + h(p, x, dom));
        end
        out = real(out);
    end

% and now we are finally ready to assemble the boundary condition matrix
% and the rhs
boundary_mat = zeros(4, 4);
for i = 1:4
   boundary_mat(1, i) = phi(3, i);
   
   boundary_mat(2, i) = phi(1, i)*L1212*r(i) - phi(3, i)*L1221*2*w_c;
   
   boundary_mat(3, i) = phi(1, i)*exp(r(i)*L1);
   boundary_mat(3, i+4) = -phi(1, i)*exp(r(i)*L1);
      
   boundary_mat(4, i) = phi(3, i)*exp(r(i)*L1);
   boundary_mat(4, i+4) = -phi(3, i)*exp(r(i)*L1);
   
   boundary_mat(5, i) = phi(1, i)*L1212*r(i)*exp(r(i)*L1) - phi(3, i)*L1221*(2*w_c)*exp(r(i)*L1);
   boundary_mat(5, i+4) = -k*(phi(1, i)*L1212*r(i)*exp(r(i)*L1) - phi(3, i)*L1221*(2*w_c)*exp(r(i)*L1));
   
   boundary_mat(6, i) = phi(1, i)*L1122*(2*w_c)*exp(r(i)*L1) + phi(3, i)*L2222*r(i)*exp(r(i)*L1);
   boundary_mat(6, i+4) = -k*(phi(1, i)*L1122*(2*w_c)*exp(r(i)*L1) + phi(3, i)*L2222*r(i)*exp(r(i)*L1));
   
   boundary_mat(7, i+4) = phi(1, i)*L1212*r(i)*exp(r(i)*L2) - phi(3, i)*L1221*(2*w_c)*exp(r(i)*L2);
   boundary_mat(8, i+4) = phi(1, i)*L1122*(2*w_c)*exp(r(i)*L2) + phi(3, i)*L2222*r(i)*exp(r(i)*L2);
end
rhs = [ 0    ;
        -F1(0, 1);
        -w1p(L1, 1);
        -w2p(L1, 1);
        -L1212*(w1p_1(L1, 1) - k*w1p_1(L1, 2)) + L1221*(2*w_c)*(w2p(L1, 1) - k*w2p(L1, 2)) - F1(L1, 1) + k*F1(L1, 2);
        -L2222*(w2p_1(L1, 1) - k*w2p_1(L1, 2)) - L1122*(2*w_c)*(w1p(L1, 1) - k*w1p(L1, 2)) - F2(L1, 1) + k*F2(L1, 2);
        -F1(L2, 2) - L1212*w1p_1(L2, 2) + L1221*(2*w_c)*w2p(L2, 2) ;
        -F2(L2, 2) - L1122*(2*w_c)*w1p(L2, 2) - L2222*w2p_1(L2, 2) ];
    
    
    
    
C = boundary_mat\rhs;
PD.C = C;
PD.phi1 = phi(1, :);
PD.phi3 = phi(3, :);
PD.phi_inv_T_2 = phi_inv(:, 2);
PD.phi_inv_T_4 = phi_inv(:, 4);
PD.r = r;

end