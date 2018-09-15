function [] = generate_piecewiseConstant_inputFiles()

k = 2:0.03:5;
fprintf('\n');
for i = 1:length(k)
   PD = do_asymptotics_pieceConst_and_print_file(0.33, 0.9,  k(i), 0.5, i); 
   if(PD.good == 0)
       fprintf('_');
   else       
       fprintf('|');
   end
end
fprintf('\n');


end