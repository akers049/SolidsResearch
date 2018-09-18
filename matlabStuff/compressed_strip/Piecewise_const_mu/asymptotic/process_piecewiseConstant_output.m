function [] = process_piecewiseConstant_output()

exp_files = dir('/home/andrew/dealii/SolidsResearch/CompressedStrip_Project/asymptoticCalculation/output/piecewiseConstant_mu_old/asymptotic*');

data_mat = [];
for i = 1:length(exp_files)
    nextData = importdata([exp_files(i).folder,'/', exp_files(i).name]);
    data_mat = [data_mat, nextData];
end

% angle = atan2d(data_mat(3, :), data_mat(2, :));
% figure(2);
% plot(data_mat(1, :), angle, '.');
% grid on;
figure(2)
yyaxis left
plot(data_mat(1, :), data_mat(2, :), '.');
ylim([-400 400])
ylabel('\lambda_2', 'rot', 1, 'Fontsize', 18)


yyaxis right 
plot(data_mat(1, :), data_mat(3, :), '.');
ylim([-1200, 1200])
ylabel('\Lambda_2', 'rot', 1, 'Fontsize', 18)

%% align zero for left and right
xlim([2.2 3.4])
% 
% yyaxis right; ylimr = get(gca,'Ylim');ratio = ylimr(1)/ylimr(2);
% yyaxis left; yliml = get(gca,'Ylim');
% if yliml(2)*ratio<yliml(1)
%     set(gca,'Ylim',[yliml(2)*ratio yliml(2)])
% else
%     set(gca,'Ylim',[yliml(1) yliml(1)/ratio])
% end

xlabel('$\frac{\mu_1}{\mu_2}$', 'Interpreter', 'latex', 'Fontsize', 18)
grid on;

end