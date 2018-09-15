function [] = process_exponential_output()

exp_files = dir('/home/andrew/dealii/SolidsResearch/CompressedStrip_Project/asymptoticCalculation/output/exponential_mu/asymptotic*');

data_mat = [];
for i = 1:length(exp_files)
    nextData = importdata([exp_files(i).folder,'/', exp_files(i).name]);
    data_mat = [data_mat, nextData];
end

angle = atan2d(data_mat(3, :), data_mat(2, :));
figure(1)
yyaxis left
ylim([-20 40])
plot(data_mat(1, :), data_mat(2, :), '.');
ylabel('\lambda_2', 'rot', 1 , 'Fontsize', 18)



yyaxis right 
plot(data_mat(1, :), data_mat(3, :), '.');
ylim([-500, 400])
ylabel('\Lambda_2', 'rot', 1,  'Fontsize', 18)
%% align zero for left and right
yyaxis right; ylimr = get(gca,'Ylim');ratio = ylimr(1)/ylimr(2);
yyaxis left; yliml = get(gca,'Ylim');
if yliml(2)*ratio<yliml(1)
    set(gca,'Ylim',[yliml(2)*ratio yliml(2)])
else
    set(gca,'Ylim',[yliml(1) yliml(1)/ratio])
end
xlim([1.5, 5.5])
xlabel('\kappa', 'Fontsize', 18);

grid on;
end