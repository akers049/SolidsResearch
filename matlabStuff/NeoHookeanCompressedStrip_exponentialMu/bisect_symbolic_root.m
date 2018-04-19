function [value] = bisect_symbolic_root(symbolicFun, a, b, tol, maxIter)
% does the bisect method cause matlab sucks. The a and b are the intial
% guessed that have to "straddle" the root. everything else is self
% explainatory.

N = 1;
while N < maxIter
    c = (a+b)/2;
    if (b - a)/2 < tol
        break;
    end
    N = N +1;
    f_a = real(double(symbolicFun(a)));
    f_c = real(double(symbolicFun(c)));
    if( sign(f_a) == sign(f_c))
        a = c;
    else
        b = c;
    end
end
value = c;
end