function [] = process_piecewiseConstant_output(indx)

exp_files = dir('/home/andrew/Research/MinnesotaStuff/SolidsResearch/CompressedStrip_Project/asymptoticCalculation/output/piecewiseConstant_mu/piecewiseConstInputFiles_dispControl/asymptotic*');
fclose('all')
data_mat = [];
for i = 1:length(exp_files)
%     nextData = importdata([exp_files(i).folder,'/', exp_files(i).name])
    next_data = cell2mat(textscan(fopen([exp_files(i).folder,'/', exp_files(i).name]),'%f\n%f\n%f\n%f'));
    data_mat = [data_mat, next_data'];
    fclose('all');
end

% angle = atan2d(data_mat(3, :), data_mat(2, :));
% figure(2);
% plot(data_mat(1, :), angle, '.');
% grid on;
% figure(2)
% yyaxis left
% plot(data_mat(1, :), data_mat(2, :), '.');
% ylim([-400 400])
% ylabel('\lambda_2', 'rot', 1, 'Fontsize', 18)
dispControlStable = data_mat(3, :) > 0;
loadControlStable = data_mat(4, :) > 0;
moduli = data_mat(1, :);
L1 = data_mat(2, :);

figure(1);
hold on;
colors  = {[247 0 53]/255, [0 217 49]/255};
for i = 1:length(dispControlStable)
    if (isnan(data_mat(3, i)) || L1(i) > 0.95)
        continue;
    end
    plot(  1 - L1(i), moduli(i), 'x', 'Color', colors{dispControlStable(i) + 1}, 'MarkerSize', 5)
end
set(0,'defaulttextinterpreter','latex');

set(groot,'defaultAxesTickLabelInterpreter','latex');

set(groot,'defaulttextinterpreter','latex');

set(groot,'defaultLegendInterpreter','latex');
ax = gca;
ax.YAxis.FontSize = 12;
ax.XAxis.FontSize = 12;
xlabel('$\frac{t}{L}$',...
    'interpreter','latex','fontsize',18)
ylabel('$\frac{\mu_f}{\mu_s}      $',...
    'Interpreter','latex','FontSize',18, 'rot', 1, 'VerticalAlignment','middle', 'HorizontalAlignment', 'right')
axis([0.0 0.53 1.9 2.9])
set(gcf, 'Color', 'None');
set(gca, 'Color', 'None');
fileName = ['disp_control_' int2str(indx), '.pdf'];
export_fig(fileName, '-pdf')

figure(2);
hold on;
for i = 1:length(loadControlStable)
    if (isnan(data_mat(4, i)))
        continue;
    end
    plot(1 - L1(i), moduli(i), 'x', 'Color', colors{loadControlStable(i) + 1}, 'MarkerSize', 5)
end
axis([0 0.53 1.5 8.5])
set(0,'defaulttextinterpreter','latex');

set(groot,'defaultAxesTickLabelInterpreter','latex');

set(groot,'defaulttextinterpreter','latex');

set(groot,'defaultLegendInterpreter','latex');

ax = gca;
ax.YAxis.FontSize = 12;
ax.XAxis.FontSize = 12;

xlabel('$\frac{t}{L}$',...
    'interpreter','latex','fontsize',18)
ylabel('$\frac{\mu_f}{\mu_s}$        ',...
    'Interpreter','latex','FontSize',18, 'rot', 1, 'VerticalAlignment','middle', 'HorizontalAlignment', 'right')
fileName_load = ['load_control_' int2str(indx), '.pdf'];
set(gcf, 'Color', 'None');
set(gca, 'Color', 'None');
export_fig(fileName_load, '-pdf')


bothStable = dispControlStable + loadControlStable;
figure(3);
hold on;
colors  = {[247 0 53]/255, [0 85 255]/255,[0 217 49]/255};
for i = 1:length(dispControlStable)
    if (isnan(data_mat(3, i)) || L1(i) > 0.95)
        continue;
    end
    plot(  1 - L1(i), moduli(i), 'x', 'Color', colors{bothStable(i) + 1}, 'MarkerSize', 5)
end
set(0,'defaulttextinterpreter','latex');

set(groot,'defaultAxesTickLabelInterpreter','latex');

set(groot,'defaulttextinterpreter','latex');

set(groot,'defaultLegendInterpreter','latex');
ax = gca;
ax.YAxis.FontSize = 12;
ax.XAxis.FontSize = 12;
xlabel('$\frac{t}{L}$',...
    'interpreter','latex','fontsize',18)
ylabel('$\frac{\mu_f}{\mu_s}      $',...
    'Interpreter','latex','FontSize',18, 'rot', 1, 'VerticalAlignment','middle', 'HorizontalAlignment', 'right')
axis([0.0 0.53 1.9 2.9])
set(gcf, 'Color', 'None');
set(gca, 'Color', 'None');
fileName = ['phase_diagram_both_' int2str(indx), '.pdf'];
export_fig(fileName, '-pdf')
% yyaxis right 
% plot(data_mat(1, :), data_mat(3, :), '.');
% ylim([-1200, 1200])
% ylabel('\Lambda_2', 'rot', 1, 'Fontsize', 18)
% 
% %% align zero for left and right
% xlim([2.2 3.4])
% % 
% % yyaxis right; ylimr = get(gca,'Ylim');ratio = ylimr(1)/ylimr(2);
% % yyaxis left; yliml = get(gca,'Ylim');
% % if yliml(2)*ratio<yliml(1)
% %     set(gca,'Ylim',[yliml(2)*ratio yliml(2)])
% % else
% %     set(gca,'Ylim',[yliml(1) yliml(1)/ratio])
% % end
% 
% xlabel('$\frac{\mu_1}{\mu_2}$', 'Interpreter', 'latex', 'Fontsize', 18)
% grid on;

end