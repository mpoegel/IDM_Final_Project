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

[num_exp, num_slopes] = size(Edata);
[num_expr, num_slopesr] = size(Rdata);

load('Experiments.mat');

%Here is the content of Names. Each one is a vector 
Abbr 
ExperimentNo
Protein 
Concentration 
Kelvin


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
for i = 1:num_exp
    cc = text(principal_coordinates(i,1),principal_coordinates(i,2),Abbr(i));
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
for i = 1:num_exp
    cc = text(principal_coordinates(i,3),principal_coordinates(i,4),Abbr(i));
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
for i = 1:num_exp
    cc = text(principal_coordinates(i,5),principal_coordinates(i,6),Abbr(i));
end


%% Clustering using K-means on Data by Experiment 

X = principal_coordinates(:,1:2);

% find the elbow in the graph
figure
title(' ');
xlabel('');
ylabel('');
K = [];
OBJ = [];
% loop over different k values
for k = 1:10
    [IDX, C] = kmeans(X, k);
    [Objective, DBI] = getDBobj(X, IDX, C);
    K(end+1,:) = k;
    OBJ(end+1,:) = Objective;
end
plot(K,OBJ,'k.-');

% graph the data using the best K

k = 3;
[IDX, C] = kmeans(X, k);

figure
grid on
title('Data by Experiment with K Clusters');
xlabel('Principal Coordinate 1');
ylabel('Principal Coordinate 2');
% Set the scale of the graph.  
s=max(max(abs(X)))*1.1;
axis([-s s -s s]);

color_sym = {'r','b','m','c','y','k','g',};

for i = 1:num_exp
    cc = text(X(i,1),X(i,2),Abbr(i));
    color_choice = color_sym{ mod(IDX(i),6)+1 };
    set(cc,'Color',color_choice);
end



%% Histograms of left and right slopes by Experiment

% we can print out all 19 experimentsa but I don't know what it means exactly and they're squished,
% maybe we need to think of another way to plot the histograms or figure
% out the slope hist code to properly edit it. The slopes generally appear
% to measure similarly. 
assigned = [8,9,11,15];


for i=1:4
    
    counter = 0;
    leftsl_ = [];
    rightsl_ = [];
    
    for j = 1:112
        if mod(counter,2) == 0
            leftsl_(:,end+1) = Edata(assigned(i),j);
    
        else
            rightsl_(:,end+1) = Edata(assigned(i),j);
       
        end
    
        counter = counter + 1;
    end
    
    output = int2str(assigned(i));
   %% blanks is suppose to be a space character but it isn't working relatively
   
    title = strcat('Experiment', blanks(10), output);
    SlopeHist( leftsl_, rightsl_, title)
end



%% PCA of Data by Residue


[eigenvectors, principal_coordinates, D] = pca(Rdata);

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

%% Plot Component 1 vs Component 2
figure
grid on
title('Residue PCA By Experiment of Components 1 and 2');
xlabel('Component 1');
ylabel('Component 2');
% Set the scale of the graph.  
s=max(max(abs(principal_coordinates1(:,1:2))))*1.1;
axis([-s s -s s]);
%The text command print 
for i = 1:num_expr
    cc = text(principal_coordinates1(i,1),principal_coordinates1(i,2),Protein(i));
end

% Plot Component 3 vs Component 4
figure
grid on
title('Residue PCA By Experiment of Components 3 and 4');
xlabel('Component 3');
ylabel('Component 3');
% Set the scale of the graph.  
s=max(max(abs(principal_coordinates1(:,3:4))))*1.1;
axis([-s s -s s]);
%The text command print 
for i = 1:num_expr
    cc = text(principal_coordinates1(i,3),principal_coordinates1(i,4),Protein(i));
end

% Plot Component 5 vs Component 6
figure
grid on
title('Residue PCA By Experiment of Components 5 and 6');
xlabel('Component 5');
ylabel('Component 6');
% Set the scale of the graph.  
s=max(max(abs(principal_coordinates1(:,5:6))))*1.1;
axis([-s s -s s]);
%The text command print 
for i = 1:num_expr
    cc = text(principal_coordinates1(i,5),principal_coordinates1(i,6),Protein(i));
end

%% Clustering using K-means on Data by Residue



%% Bar Graph of for left and right slopes of a residue





