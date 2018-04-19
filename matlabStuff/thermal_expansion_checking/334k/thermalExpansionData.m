clear
clc

%% Thermal Expansion Data

K = 99.4973e9;  % the fitted value from the paper



rawPressures = [-53.398986;
                579.17102 ;
                1227.6955 ;
                1892.6081 ;
                2574.3564;
                3273.4024]*1e5;  % in Pascals
deltaT = [0:10:50]';
deltaP = rawPressures - rawPressures(1);
deltaP_by_3K = deltaP/(3*K);

lin_alpha = lsqnonneg(deltaT(1:end-1), deltaP_by_3K(1:end-1));

% visualize data
figure(1)

t = 0:0.01:40;
y = lin_alpha*t;
plot(t, y, 'k', 'Linewidth', 2);
hold on
plot(deltaT(1:end-1), deltaP_by_3K(1:end-1), 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
grid on

xlabel('\DeltaT (K)');
ylab = ylabel('$\frac{\Delta P}{3\mathcal{K}}$');
set(ylab,'Interpreter','latex','FontSize', 18);

leg = legend('Linear fit with \alpha = 21.825e-6 K^{-1}', 'Values from simulation', 'Location', 'NorthWest') ;
set(leg, 'FontSize', 12);
title('\DeltaP vs \DeltaT for simulated and linear fit');