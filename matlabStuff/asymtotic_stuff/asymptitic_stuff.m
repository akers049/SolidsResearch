clear
clc
close all

syms kappa x1 x2 nu l1 l2 Ai ai w AiBi;

syms A1 A1B1 A2 A2B2 a1 a2;

u1 = -Ai*exp(ai*x2)*sin(w*x1);
u2 = AiBi*exp(ai*x2)*cos(w*x1);

u1_1 = diff(u1, x1);
u1_2 = diff(u1, x2);

u2_1 = diff(u2, x1);
u2_2 = diff(u2, x2);

syms M111111 M222222 M111122 M112222 M111212 M112121 ...
     M112112 M12112 M121121 M122212 M122221 M221122 ...
     M221111 M221212 M222121 M221221 M211112 M211121 M212221 M212212;
 
 % for i = 1
 term_1_1 = M111111*(u1_1)^2 + 2*M111122*u1_1*u2_2 + M112222*(u2_2)^2 + ...
     M111212*(u1_2)^2 + M112121*(u2_1)^2 + 2*M112112*u1_2*u2_1;
 term_1_1 = simplify(term_1_1);

 term_1_2 = 2*(M12112*u1_1*u1_2 + M121121*u1_1*u2_1 + ...
     M122212*u2_2*u1_2 + M122221*u2_2*u2_1);
 term_1_2 = simplify(term_1_2);
 
 rhs_1 = -( diff(term_1_1, x1) + diff(term_1_2, x2) + kappa*term_1_2);
 rhs_1 = simplify(rhs_1)
 
 % for i = 2
 term_2_1 = 2*(M211112*u1_1*u1_2 + M211121*u1_1*u2_1 + ...
     M212221*u2_2*u2_1 + M212212*u2_2*u1_2);
 term_2_1 = simplify(term_2_1);
 
 term_2_2 = M222222*(u2_2)^2 + 2*M221122*u1_1*u2_2 + M221111*(u1_1)^2 + ...
     M221212*(u1_2)^2 + M222121*(u2_1)^2 + 2*M221221*u1_2*u2_1;
 term_2_2 = simplify(term_2_2);
 
 rhs_2 = -(diff(term_2_1, x1) + diff(term_2_2, x2) + kappa*term_2_2);
 rhs_2 = simplify(rhs_2)
 
