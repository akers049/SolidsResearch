clear
clc
close all

syms a b c d e1 e2 e3 ee1 ee2 ee3 f g h k l;

syms L1111 L1212 L1122 L2112 L2222 dub E1 E1_1 E1_11 E2 E2_1 E2_11 ep;

a = L1212;
c  = -L1111*dub*dub;
d = -dub*(L2112 + L1122);
e1 = -ep*E1;
e2 = -ep*E1_1;
e3 = -ep*E1_11;

f = dub*(L1122+L2112);
h = L2222;
l = -dub*dub*L1212;
ee1 = -ep*E2;
ee2 = -ep*E2_1;
ee3 = -ep*E2_11;



systemMat = [ 0 0 a 0 c 0 0 0 d 0 e1; 
              0 a 0 c 0 0 0 d 0 0 e2;
              a 0 c 0 0 0 d 0 0 0 e3;
              0 0 0 f 0 0 0 h 0 l ee1; 
              0 0 f 0 0 0 h 0 l 0 ee2;
              0 f 0 0 0 h 0 l 0 0 ee3];
          

systemMat(3, :) = systemMat(3, :) - systemMat(3, 7)/ systemMat(5, 7)*systemMat(5, :);
systemMat(3, :) = systemMat(3, :) - systemMat(3, 9)/ systemMat(1, 9)*systemMat(1, :);

systemMat(6, :) = systemMat(6, :) - systemMat(6, 2)/ systemMat(2, 2)*systemMat(2, :);
systemMat(6, :) = systemMat(6, :) - systemMat(6, 4)/ systemMat(4, 4)*systemMat(4, :);      
          
systemMat(3, :) = systemMat(3, :)*L2222/L1212;

eq1 = simplify(systemMat(3, [1:5, 11]))
eq2 = simplify(systemMat(6, [6:11]))