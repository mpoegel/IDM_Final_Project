%% Final Project
% Analysis of Protein Folding
% Matthew Poegel and Jessie Sodolo
% May 7, 2014

%% Read the data 
R=csvread('data.csv'); % this is the data organized to study by residue 
E=csvread('dataT.csv'); %this is the data organized to study by experiments

% pull off the top and bottom 
Rdata=R(2:end,2:end); % data organized by residue 
RColLabels=R(1,2:end); % Proteins corresponding to each col of R 
RRowLabels=R(2:end,1); % Residues corresponding to each row of R 
Edata=E(2:end,2:end); % Data organized by experiments 
EColLabels=E(1,2:end); % Residues corresponding to each col of E 
ERowLabels=E(2:end,1); % Experiments corresponding to each row of E 
Names=dataset('file', 'eptTable.csv'); % Get the details on the exeperiments.

%Here is the content of Names. Each one is an array 
Names.Abbr 
Names.ExperimentNo
Names.Protein 
Names.Concentration 
Names.Kelvin
