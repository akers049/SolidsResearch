function [] = generate_exponential_inputFiles()

k = 1.80130:0.000001:1.80133;
fprintf('\n');
for i = 1:length(k)
   PD = do_asymptotics_exponential_and_print_file(0.33, k(i), 0.5, i); 
   if(PD.good == 0)
       fprintf('_');
   else       
       fprintf('|');
   end
end
fprintf('\n');


end