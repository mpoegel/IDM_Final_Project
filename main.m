%% Final Project
% Analysis of Protein Folding
% Matthew Poegel and Jessie Sodolo
% May 7, 2014

close all;
clear;

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
% Names=readtable('eptTable.csv'); % Get the details on the exeperiments.

%Here is the content of Names. Each one is an array 
% Names.Abbr 
% Names.ExperimentNo
% Names.Protein 
% Names.Concentration 
% Names.Kelvin

NamesAbbr = { 'D-1.5-293'
            'D-1.8-293'
            'D-2-293'
            'D-2-288'
            'D-1.8-298'
            'D-1.8-303'
            'I-0.5-293'
            'I-0.65-293'
            'I-0.85-293'
            'I-0.5-288'
            'I-0.5-298'
            'I-0.5-303'
            'L-0.6-293'
            'L-0.75-293'
            'L-0.85-293'
            'L-0.75-288'
            'L-0.75-298'
            'L-0.75-303'
            'S-5.5-293' };

display(NamesAbbr)


%% PCA of Data of Experiment

[eigenvectors, principal_coordinates, D] = pca(Edata);

% Look at the explained variance of the eigenvalues
var_explained  = cumsum(D)/sum(D);
figure
subplot(2,1,1)
bar(D);
title('Eigenvalues')
subplot(2,1,2)
bar(var_explained)
set(gca,'YTick',0:0.1:1)
set(gca,'XTick',0:5:69)
title('Variance Explained')

% Plot Component 1 vs Component 2
figure
grid on
title('PCA By Experiment of Components 1 and 2');
xlabel('Component 1');
ylabel('Component 2');
% Set the scale of the graph.  
s=max(max(abs(principal_coordinates(:,1:2))))*1.1;
axis([-s s -s s]);
%The text command print 
for i = 1:18
    cc = text(principal_coordinates(i,1),principal_coordinates(i,2),NamesAbbr(i));
end

% Plot Component 3 vs Component 4
figure
grid on
title('PCA By Experiment of Components 3 and 4');
xlabel('Component 3');
ylabel('Component 4');
% Set the scale of the graph.  
s=max(max(abs(principal_coordinates(:,3:4))))*1.1;
axis([-s s -s s]);
%The text command print 
for i = 1:18
    cc = text(principal_coordinates(i,3),principal_coordinates(i,4),NamesAbbr(i));
end

% Plot Component 5 vs Component 6
figure
grid on
title('PCA By Experiment of Components 5 and 6');
xlabel('Component 5');
ylabel('Component 6');
% Set the scale of the graph.  
s=max(max(abs(principal_coordinates(:,5:6))))*1.1;
axis([-s s -s s]);
%The text command print 
for i = 1:18
    cc = text(principal_coordinates(i,5),principal_coordinates(i,6),NamesAbbr(i));
end


%% Clustering using K-means on Data by Experiment 



%% Histograms of left and right slopes by Experiment



%% PCA of Data by Residue



%% Clustering using K-means on Data by Residue



%% Bar Graph of for left and right slopes of a residue





