function [d3W_dFdFdF] = get_d3Wc_dFdFdF(mu0, nu, lambda)

% give lambda1, the stretch in the x1 direction
lambda1 = 1 - lambda;

% solve for lambda2, the stretch in the x2 direction (use the positive
% solution).
lambda2_sols = roots([(mu0 + lambda1^2*2*nu*mu0/(1 - nu)), -(2*mu0*nu/(1 - nu))*lambda1, -mu0]);
lambda2 = (lambda2_sols(find(double(lambda2_sols >= 0))));

F_inv = [1/lambda1 0; 0 1/lambda2];
II_F = lambda1*lambda2;
II_C = II_F*II_F;

d3W_dFdFdF = zeros(2,2,2,2,2,2);
for i = 1:2
    for j = 1:2
        for k = 1:2
            for l = 1:2
                for m = 1:2
                    for n = 1:2
                        d3W_dFdFdF(i,j,k,l,m,n) = (-1 + (2*nu/(1-nu))*(II_C - II_F))*...
                                (F_inv(j, m)*F_inv(n, k)*F_inv(l, i) + F_inv(j, k)*F_inv(l, m)*F_inv(n, i)) + ...
                            (-4*nu/(1-nu))*(II_C - 0.5*II_F)*...
                                (F_inv(n, m)*F_inv(j, k)*F_inv(l, i) + F_inv(l, m)*F_inv(n, k)*F_inv(j, i) + F_inv(l, k)*F_inv(j, m)*F_inv(n, i)) + ...
                            (8*nu/(1-nu))*(II_C - 0.25*II_F)*F_inv(n, m)*F_inv(l, k)*F_inv(j, i);
                    end
                end
            end
        end
    end
end

d3W_dFdFdF = mu0*d3W_dFdFdF;
end