clear
clc
close all
data_cell = cell(0,0);
files = dir('../CompressedStrip_Project/output/1_6_1/bloch_eigenvalues/bloch_eigenvalues-*');
for i = 1:length(files)
   nextFilePath = [files(i).folder , '/' , files(i).name];
   %index = str2double(files(i).name(13:end));
   
   nextDataVect = importdata(nextFilePath);
   data_cell{i} = nextDataVect; 

    
end
save('eigenvalues_path_data.mat', 'data_cell');

%% 
load('eigenvalues_path_data.mat');
figure(1);
hold on;
for i=1:length(data_cell)
    for k = 1:length(data_cell{i})
        if  data_cell{i}(k) < 0.0
           plot(i, data_cell{i}(k), 'r.', 'MarkerSize', 7); 
        else
           plot(i, data_cell{i}(k), 'k.', 'MarkerSize', 7); 
        end
    end

end
grid on;