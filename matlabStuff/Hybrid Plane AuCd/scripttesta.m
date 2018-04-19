clear
clc
for i = 0:100000
C = [1 2 12; 44 23 9; 1 1 9];
[a, b] = calcAndOrderEigs(C);
end