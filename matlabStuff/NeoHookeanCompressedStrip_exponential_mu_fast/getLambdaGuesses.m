clear
clc
close all
load('1_6_frequencies_and_lambdas.mat');
min_index = find(min(solutionVector) == solutionVector);
min_wavelength = 2*pi/w_vector(min_index);

wavelength_vector = fliplr(2*pi*w_vector.^(-1));
solutionVector = fliplr(solutionVector);

primitiveCellLength = min_wavelength/30;
fid = fopen('lambda_guesses_30.dat', 'w');

for i = 9:3:300
    next_lambda = interp1(wavelength_vector, solutionVector, i*primitiveCellLength);
    fprintf(fid, '%g     %f\n', i, next_lambda);
end

fclose(fid);
%% 

load('1_6_frequencies_and_lambdas.mat');

min_index = find(min(solutionVector) == solutionVector);
min_wavelength = 2*pi/w_vector(min_index);
primitiveCellLength = min_wavelength/30;

wavelength_vector = fliplr(2*pi./w_vector);
solutionVector = fliplr(solutionVector);




lambda_crits = load('lambda_crits_30.dat');
figure(1);
plot(wavelength_vector, solutionVector, 'lineWidth', 1.3);
hold on;
plot(lambda_crits(:, 1)*primitiveCellLength, lambda_crits(:, 2), 'x');
axis([0 20 0.3 0.8]);
grid on;
xlabel('Wavelength')
ylabel('\lambda_c     ', 'rot', 1)

