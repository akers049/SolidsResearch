function [X, Y, Lambda_c, Wavelength_c] = generate_piecewiseConstant_inputFiles()

tic
% L1 = 0.5:0.01:0.99;

% L1 = logspace(log10(0.03), log10(0.5), 50);

L1 = logspace(log10(0.02), log10(0.5), 80);
L1 = 1.0 - L1;

k = logspace(log10(1.01), log10(50), 80); %2.0:0.3333:8; % 2.0:0.05:2.8;  %1.5:0.05:2.5;
klen = length(k);
Llen = length(L1);
[X, Y] = meshgrid(L1, k);
Lambda_c = zeros(size(X));
Wavelength_c = Lambda_c;
fprintf('\n');
PDS = cell(size(X));
parfor j = 1:Llen
    size(k);
    for i = 1:klen

        
        fprintf('start %u %u\n', i, j);
%         PDS{i,j} = struct;
%         PDS{i,j}.lambda_c = 1.0;
%         PDS{i,j}.wavelength_c = 2.0;
%         PDS{i,j}.good = 1;
       PDS{i, j} = do_asymptotics_pieceConst_and_print_file(0.33, L1(j),  k(i), 0.1, i*Llen + j, true); 
       Lambda_c(i, j) = PDS{i,j}.lambda_c;
       Wavelength_c(i,j) = PDS{i, j}.wavelength_c;
       if(PDS{i,j}.good == 0)
           fprintf('_');
       else       
           fprintf('%u %u\n', i, j);
       end
    end
end
filename= 'criticalInfor_over_param_space_range_good.mat';
save(filename, 'X', 'Y', 'Lambda_c', 'Wavelength_c');
fprintf('\n');
toc

end