function [] = plot_bifurcation_diagram(folder_prefix)

folders = dir(['/home/andrew/dealii/SolidsResearch/CompressedStrip_Project/output/', folder_prefix,'*']);

load_data = cell(0,0);

for i = 1:length(folders)
    
    if (strcmp(folders(i).name((end-1):end), '_1'))
        load_mat = import_load_info([folders(i).folder, '/', folders(i).name ,'/load_info/load_info1.txt']);
        load_data{i, 1} = load_mat(1:600, 1);
        load_data{i, 2} = load_mat(1:600, 3);
        load_data{i, 3} = load_mat(1:600, 4);
        
        
        load_mat = import_load_info([folders(i).folder, '/', folders(i).name ,'/load_info/load_info4.txt']);
        load_data{length(folders)+1, 1} = load_mat(:, 1);
        load_data{length(folders)+1, 2} = load_mat(:, 3);
        load_data{length(folders)+1, 3} = load_mat(:, 4);
        
    else
        load_mat = import_load_info([folders(i).folder, '/', folders(i).name ,'/load_info/load_info2.txt']);
        load_data{i, 1} = load_mat(:, 1);
        load_data{i, 2} = load_mat(:, 3);
        load_data{i, 3} = load_mat(:, 4);
    end

end

figure(1);
grid on
hold on

for i = 1:(length(folders) + 1)
   plot(load_data{i, 1}, load_data{i, 2}, 'LineWidth', 1.2);  
end
title('Bifurcation Diagram for Piecewise-constant $\mu$', 'interpreter', 'latex')
xlabel('\lambda');
ylabel('\Lambda      ', 'rot', 1);


figure(2)
grid on
hold on
for i = 1:(length(folders) + 1)
   plot(load_data{i, 1}, load_data{i, 3}, 'LineWidth', 1.2);  
end
ylabel('|u|_2     ', 'rot', 1);
xlabel('\lambda');
title('Bifurcations with Displacement Magnitude', 'interpreter', 'latex')

end