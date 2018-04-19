clear
clc
close all

lambda = 0.1:0.01:0.95;
cows = zeros(size(lambda));
for i = 1:length(cows)
   cows(i) = evaluate_determinant_M(0.1, 1.7, lambda(i));
    
end

plot(lambda, cows)
axis([0.1, 0.95, -10, 10])