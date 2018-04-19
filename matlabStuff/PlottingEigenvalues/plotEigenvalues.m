function [] = plotEigenvalues(stepsPerOutput)

cd output\;
files = dir('eigenvalues-*');
figure(1);
hold on;
for i = 1:length(files);
   nextFileName = files(i).name;
   A = importdata(nextFileName);
   
   eigenValues = A.data;
   stepNum = stepsPerOutput*str2double(nextFileName(13:(end)));
   
   plot(stepNum, eigenValues, 'k.', 'MarkerSize', 5);
end
grid on;



end