clear
clc
close all

syms a b c d e f g h k l;

syms L1111 L1212 L1122 L2112 L2222 eta1 eta3 m w;

a = L1212;
b = m*L1212;
c  = -w*w*L1111;
d = -w*(L2112 + L1122);
e = -w*m*L2112;
f = w*(L1122+L2112);
g = w*m*L1122;
h = L2222;
k = L2222*m;
l = -w*w*L1212;


systemMat = [ 0 0 a b c 0 0 0 d e; 
              0 a b c 0 0 0 d e 0;
              a b c 0 0 0 d e 0 0;
              0 0 0 f g 0 0 h k l; 
              0 0 f g 0 0 h k l 0;
              0 f g 0 0 h k l 0 0];
          
systemMat(4, :) = systemMat(4, :) - systemMat(4, 8)/ systemMat(2, 8)*systemMat(2, :)
systemMat(4, :) = systemMat(4, :) - systemMat(4, 9)/ systemMat(1, 9)*systemMat(1, :)
systemMat(3, :) = systemMat(3, :) - systemMat(3, 7)/ systemMat(5, 7)*systemMat(5, :)
systemMat(3, :) = systemMat(3, :) - systemMat(3, 8)/ systemMat(2, 8)*systemMat(2, :)
systemMat(3, :) = systemMat(3, :) - systemMat(3, 9)/ systemMat(1, 9)*systemMat(1, :)
systemMat(3, :) = systemMat(3, :) - systemMat(3, 10)/ systemMat(4, 10)*systemMat(4, :)

systemMat(4, :) = [0 0 0 f g 0 0 h k l];
y1System = systemMat(3, 1:5);
y1System = simplify(y1System)

systemMat(1, :) = systemMat(1, :) - systemMat(1, 3)/ systemMat(5, 3)*systemMat(5, :)
systemMat(1, :) = systemMat(1, :) - systemMat(1, 4)/ systemMat(4, 4)*systemMat(4, :)
systemMat(6, :) = systemMat(6, :) - systemMat(6, 2)/ systemMat(2, 2)*systemMat(2, :)
systemMat(6, :) = systemMat(6, :) - systemMat(6, 3)/ systemMat(5, 3)*systemMat(5, :)
systemMat(6, :) = systemMat(6, :) - systemMat(6, 4)/ systemMat(4, 4)*systemMat(4, :)
systemMat(6, :) = systemMat(6, :) - systemMat(6, 5)/ systemMat(1, 5)*systemMat(1, :)

y2System = systemMat(6, 6:end);
y2System = simplify(y2System)