clear
clc
close all
%% 
 % directory_name = '~/dealii/SolidsResearch/CompressedStrip_Project/output/nonUniform_3_0_1';
%  directory_name = '~/dealii/SolidsResearch/CompressedStrip_Project/output/coarse_pieceConst_10_0_9_1';
  directory_name = '~/dealii/SolidsResearch/CompressedStrip_Project/output/nonUniform__6_0_1';

indx = [0];


data_cell = cell(0,0);
files = dir([directory_name, '/bloch_eigenvalues/bloch_eigenvalues', int2str(indx),  '-*']);
for i = 1:length(files)
   nextFilePath = [files(i).folder , '/' , files(i).name];
   nextDataVect = importdata(nextFilePath);
   
   
   
   if (length(data_cell) < nextDataVect(1))
      data_cell(nextDataVect(1)) = {cell(0,0)}; 
   end
   
   data_cell{nextDataVect(1)} = [data_cell{nextDataVect(1)}, {nextDataVect(2:end)}]; 

    
end
save('eigenvalues_data.mat', 'data_cell');



%%

files = dir([directory_name, '/load_info/load_info0*']);

nextFilePath = [files(1).folder , '/' , files(1).name];
nextDataVect = importdata(nextFilePath);
load_data = nextDataVect.data;
save('load_data.mat', 'load_data');

%%
% close all;
load('eigenvalues_data.mat');
max = length(data_cell);
f = figure(1);

ax2 = axes('Parent',f,'position',[0.13 0.25  0.77 0.27]);
ax = axes('Parent',f,'position',[0.13 0.60  0.77 0.27]);

f.UserData = [{ax}, {ax2}];


b = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],...
              'value',1.0, 'min',1, 'max',max);
bgcolor = f.Color;
bl1 = uicontrol('Parent',f,'Style','text','Position',[50,54,23,23],...
                'String','1','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',f,'Style','text','Position',[500,54,23,23],...
                'String',string(max),'BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',f,'Style','text','Position',[240,25,100,23],...
                'String','step','BackgroundColor',bgcolor);
b.Callback = @(es,ed) slider2_Callback(es, ed); 
            

            