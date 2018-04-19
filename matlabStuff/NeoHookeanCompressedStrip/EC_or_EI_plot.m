clear
clc
close all


lambdas = 0.01:0.0001:0.9;
EI_or_EC = zeros(size(lambdas));
for i =1:length(lambdas)
    EI_or_EC(i) = determine_EC_or_Ei(lambdas(i), 0.33);
end
plot(lambdas, EI_or_EC);