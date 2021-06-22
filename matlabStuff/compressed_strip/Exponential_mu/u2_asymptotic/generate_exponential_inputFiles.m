function [] = generate_exponential_inputFiles()
tic
k = 0.5:0.02:8.5;
fprintf('\n');
for i = 1:length(k)
   PD = do_asymptotics_exponential_and_print_file(0.33, k(i), 0.1, i); 
   if(PD.good == 0)
       fprintf('_');
   else       
       fprintf('|');
   end
end
fprintf('\n');
toc

end