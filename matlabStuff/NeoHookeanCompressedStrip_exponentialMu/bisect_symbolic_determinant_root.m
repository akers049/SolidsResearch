function [value] = bisect_symbolic_determinant_root(symbolicMatrixFunction, a, b, tol, maxIter)
% does the bisect method cause matlab sucks. The a and b are the intial
% guessed that have to "straddle" the root. The symbolicMatrixFunction
% is a symbolic function that takes in one argument (lambda) and returns the 
% symbolic matrix that the root wants to be found for. everything else is self
% explainatory.

N = 1;
while N < maxIter
    c = (a+b)/2;
    if (b - a)/2 < tol
        break;
    end
    N = N +1;
    f_a = real(det(double(symbolicMatrixFunction(a))));
    f_c = real(det(double(symbolicMatrixFunction(c))));
    if( sign(f_a) == sign(f_c))
        a = c;
    else
        b = c;
    end
end
value = c;
end