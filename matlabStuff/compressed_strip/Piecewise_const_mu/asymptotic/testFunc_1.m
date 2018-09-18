function [val] = testFunc_1(PD, x1, x2)

w_c = PD.w_c;
A = PD.A;
B = PD.B;
a = PD.alphas;

v1 = 0;
v2 = 0;
t = 0;
if (x2 > PD.L1)
    t = 4;
end
for i = 1:4
   v1 = v1 + A(i + t)*exp(a(i)*x2); 
   v2 = v2 + B(i)*A(i + t)*exp(a(i)*x2); 
end

val = (v1.*v1).*(sin(w_c*x1).^2) + (v2.*v2).*(cos(w_c*x1).^2);

end