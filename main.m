%% Introduction to Data Math -- Final Project
% Analysis of Protein Folding
% Matthew Poegel and Jessie Sodolo
% May 9, 2014

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
[num_residues, num_exp2] = size(Rdata);

load('Experiments.mat');

%Here is the content of Names. Each one is a vector 
Abbr 
ExperimentNo
Protein 
Concentration 
Kelvin


%% PCA of Data of Experiment
% Using the data as organized by experiment, the analysis begins with a
% principal coordinate analysis of the experiments. We look at the
% explained variance and the plots of PC1 v PC2, PC3 v PC4, and PC5 v PC6.

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
% Analysis of the clustering the experiments using k-means. Examing k
% values of 1-10, a plot of objective v K is made to obtain the elbow. Then
% the clusters are indentified using the K value at the elbow.

% find the elbow in the graph
figure;
K = [];
OBJ = [];
% loop over different k values
for k = 1:10
    [IDX, C] = kmeans(Edata, k);
    [Objective, DBI] = getDBobj(Edata, IDX, C);
    K(end+1,:) = k;
    OBJ(end+1,:) = Objective;
end
hold on;
plot(K,OBJ,'k.-');
title('Objective vs. K Clusters');
xlabel('K Value');
ylabel('Objective');

% chart the data using the best K
k = 3;
[IDX, C] = kmeans(Edata, k);
[obj, DBI] = getDBobj(Edata, IDX, C);
% mark the elbow of the graph with a red plus
plot(k,obj,'r+');
hold off;
% count the number in each cluster
count = zeros(k,1);
for i = 1:num_exp
    count(IDX(i)) = count(IDX(i)) + 1;
end

% make a chart showing the total number in each cluster
figure
row_names = {'Cluster 1', 'Cluster 2', 'Cluster 3'};
col_names = {'Data in Cluster'};
uitable('ColumnName',col_names, 'RowName', row_names, 'Data',count, 'Position',[20 300 360 100]);


%% Histograms of left and right slopes by Experiment
% Create histograms of the slopes for experiments 1-9, 11, and 15.

%All the experiments assigned to us our in this array for usuage.
assigned = [7,8,9,11,15];

for i=1:5    
    counter = 0;
    %Separate matrix to keep the left and right slopes. 
    leftsl_ = [];
    rightsl_ = [];
    
    for j = 1:112        
        %The residues go left to right, if the place in Edata is even go
        %left if odd right.        
        if mod(counter,2) == 0
            leftsl_(:,end+1) = Edata(assigned(i),j);
    
        else
            rightsl_(:,end+1) = Edata(assigned(i),j);       
        end
    
        counter = counter + 1;
    end    
    
    % create the histogram using the SlopeHist function
    output = int2str(assigned(i));   
    T = strcat('Experiment-', output);
    SlopeHist( leftsl_, rightsl_, T);
    
end



%% PCA of Data by Residue
% Principal coordinate analysis of the residues. Examine the variance as
% explained by the eigenvalues and plot the PC1 v PC2, PC3 v PC4, and PC5 v
% PC6.

[eigenvectors, principal_coordinates1, D] = pca(Rdata);

% Look at the explained variance of the eigenvalues
var_explained  = cumsum(D)/sum(D);
figure
subplot(2,1,1)
bar(D);
title('Eigenvalues');
subplot(2,1,2)
bar(var_explained)
set(gca,'YTick',0:0.1:1)
set(gca,'XTick',0:5:69)
title('Variance Explained')

% Plot Component 1 vs Component 2
figure
grid on
title('PCA By Residue of Components 1 and 2');
xlabel('Component 1');
ylabel('Component 2');
% Set the scale of the graph.  
s=max(max(abs(principal_coordinates1(:,1:2))))*1.1;
axis([-s s -s s]);
%The text command print 
for i = 1:num_residues
    cc = text(principal_coordinates1(i,1),principal_coordinates1(i,2),int2str(RRowLabels(i)));
end

% Plot Component 3 vs Component 4
figure
grid on
title('PCA By Residue of Components 3 and 4');
xlabel('Component 3');
ylabel('Component 3');
% Set the scale of the graph.  
s=max(max(abs(principal_coordinates1(:,3:4))))*1.1;
axis([-s s -s s]);
%The text command print 
for i = 1:num_residues
    cc = text(principal_coordinates1(i,3),principal_coordinates1(i,4),int2str(RRowLabels(i)));
end

% Plot Component 5 vs Component 6
figure
grid on
title('PCA By Residue of Components 5 and 6');
xlabel('Component 5');
ylabel('Component 6');
% Set the scale of the graph.  
s=max(max(abs(principal_coordinates1(:,5:6))))*1.1;
axis([-s s -s s]);
%The text command print 
for i = 1:num_residues
    cc = text(principal_coordinates1(i,5),principal_coordinates1(i,6),int2str(RRowLabels(i)));
end

%% Clustering using K-means on Data by Residue
% Examine the clustering of the residue data using kmeans. Create a graph
% of objective v k value, using k values of 1-10. Find the elbow and
% examine the clusters created with this k.

% find the elbow in the graph
figure;
K = [];
OBJ = [];
% loop over different k values
for k = 1:10
    [IDX, C] = kmeans(Rdata, k);
    [Objective, DBI] = getDBobj(Rdata, IDX, C);
    K(end+1,:) = k;
    OBJ(end+1,:) = Objective;
end
hold on;
plot(K,OBJ,'k.-');
title('Objective vs. K Clusters');
xlabel('K Value');
ylabel('Objective');


% chart the data using the best K
k = 3;
[IDX, C] = kmeans(Rdata, k);
[obj, DBI] = getDBobj(Rdata, IDX, C);
% mark the elbow of the graph with a red plus
plot(k,obj,'r+');
hold off;
% count the number in each cluster
count = zeros(k,1);
for i = 1:num_residues
    count(IDX(i)) = count(IDX(i)) + 1;
end

% make a chart showing the number in each cluster
figure
row_names = {'Cluster 1', 'Cluster 2', 'Cluster 3'};
col_names = {'Data in Cluster'};
t = uitable('ColumnName',col_names, 'RowName', row_names, 'Data',count, 'Position',[20 300 360 100]);

% show the cluster that each residue is in
clusters = [RRowLabels IDX];
display(clusters)


%% Bar Graph of for left and right slopes of a residue
% create a bar graph of the left and right slopes of each residue in
% experiments 7-9, 11, and 15.

%All the experiments assigned to us our in this array for usuage.
assigned = [7,8,9,11,15];

for i = assigned;
    %loop through each experiment and store it's left slopes in one matrix and the right in another
    leftsl_ = [];
    rightsl_ = [];
    
    
    %variables used to fill in the zeros for missing residue numbers
    place = 1;
    %variable is used to obtain every residue with Edata
    j = 1;
    
    for k = 1:140       
     %Condition checks if spot which is an integer is an actual residue
     %number we have data for. If it is then we need to calculate its
     %slopes, if its not then we fill a zero in for that residue number. 
     
        if k < RRowLabels(place)
            leftsl_(:,end+1) = 0;
            rightsl_(:,end+1) = 0;
        end
        
        if k == RRowLabels(place)  
            leftsl_(:,end+1) = Edata(i,j);
            rightsl_(:,end+1) = Edata(i,j+1);

            place = place + 1;
            j = j+2;
        end
    end

    %Output the contents in each matrix.
    figure
    axis([0 149 -0.02 .015])
    output = int2str(i);   
    T = strcat('Experiment-', output);
    title(T);
    xlabel('Residue')
    ylabel('Slopes')
    hold on
    bar(leftsl_,'r');
    bar(rightsl_,'b');
    
end


%% Analysis of Slope Ratios by Experiment

% We further explored the idea of looking at the residual slopes.
% To get a better understanding of the chervon structure we came up with an
% idea. For each protein we looked up the experiments in respect it
% contained, we then summed the left side and right side separately of each
% residue through all experiments in that protein. After the summation a
% ratio of left slope divided by right slope was calculated. The ratio was
% the plotted for each residue. This process was repeated for remaining
% proteins. 


% variable used to print the left and right slopes

place = 1;

%Matrices used to store the ratio amounts calculated. 
entry = [];
secpro = [];
thridpro = [];

% variable used to fill in the zeros for missing residue numbers
spot =1;
check =1;

while place < 112 ,
    %Variable to store the residues left total for each protein
    
    % left and right for I92A
    lefttotal = 0;
    righttotal = 0;
    
    % left and right for L125A
    slefttotal = 0;
    srighttotal = 0;
   
    % left and right for D+PHS 
    lt = 0;
    rt = 0;
    
    %Variable to store each ratio protein I92A, L125A, D+PHS respectively    
    ratio = 0;
    ratio2 = 0;
    ratio3 = 0;
    
    
    %Spot is a residue number we don't have 
    if spot < RRowLabels(check); 
    entry(:,end+1) = 0; 
    secpro(:,end+1) = 0;
    thridpro(:,end+1) = 0;

    end
    
    %Condition checks if spot which is an integer is an actual residue
    %number we have data for. If it is then we need to calculate its
    %slopes.     
    if spot == RRowLabels(check); 

        % left slope for D+PHS 
        lt = Edata(1, place);
        lt = Edata(2, place) + lt;
        lt =Edata(3, place) + lt;
        lt =Edata(4, place) + lt; 
        lt =Edata(5, place) + lt;
        lt =Edata(6, place) + lt;
        lt = (lt ./ 6);

        % left slope for I92A 
        lefttotal = Edata(7, place);
        lefttotal = Edata(8, place) + lefttotal;
        lefttotal =Edata(9, place) + lefttotal;
        lefttotal =Edata(10, place) + lefttotal; 
        lefttotal =Edata(11, place) + lefttotal;
        lefttotal =Edata(12, place) + lefttotal;
        lefttotal = (lefttotal ./ 6);

        % left slope for L125A
        slefttotal = Edata(13, place);
        slefttotal = Edata(14, place) + slefttotal;
        slefttotal = Edata(15, place) + slefttotal;
        slefttotal = Edata(16, place) + slefttotal; 
        slefttotal = Edata(17, place) + slefttotal;
        slefttotal = Edata(18, place) + slefttotal;
        slefttotal = (slefttotal ./6);

        % right slope for D+PHS 
        rt = Edata(1, place+1);
        rt = Edata(2, place+1) + rt;
        rt =Edata(3, place+1) + rt;
        rt =Edata(4, place+1) + rt; 
        rt =Edata(5, place+1) + rt;
        rt =Edata(6, place+1) + rt;
        rt = (rt ./ 6);

        % right slope for I92A 
        righttotal = Edata(7, place+1);
        righttotal = Edata(8, place + 1) + righttotal;
        righttotal =Edata(9, place + 1) + righttotal;
        righttotal =Edata(10, place + 1) + righttotal; 
        righttotal =Edata(11, place + 1) + righttotal;
        righttotal =Edata(12, place + 1) + righttotal;
        righttotal = (righttotal ./6);

        % right slope for L125A 
        srighttotal = Edata(13, place + 1);
        srighttotal = Edata(14, place + 1) + srighttotal;
        srighttotal =Edata(15, place + 1) + srighttotal;
        srighttotal =Edata(16, place + 1) + srighttotal; 
        srighttotal =Edata(17, place + 1) + srighttotal;
        srighttotal =Edata(18, place + 1) + srighttotal;
        srighttotal = (srighttotal ./6);

        % Ratios calculated
        ratio = lefttotal ./ righttotal;
        ratio2 = slefttotal ./ srighttotal;
        ratio3 = lt ./ rt;

        % ratios entered into respective matrix
        thridpro(:,end+1) = ratio3;
        secpro(:,end+1) = ratio2;
        entry(:,end+1) = ratio;

        %Increment the variables to move through the residues and residue
        %number. 
        place = place + 2;
        check = check +1;
        
    end       
    
    spot = spot+1;
    
end


%Each stored matrix containning protein ratios is printed out separately
%with the corresponding protein name labled as the title.

pro3 = Protein(1);
figure
grid on
axis([-1 149 -10 20])
T = strcat('Slope ratio of Protein- ', pro3);
title(T);
xlabel('Residues')
ylabel('Slopes ratio left: right')
hold on
bar(thridpro,'g');
hold off  
       
        
    
pro = Protein(7);
figure
grid on
axis([-1 149 -10 20])
T = strcat('Slope Ratio of Protein- ', pro);
title(T);
xlabel('Residues')
ylabel('Slopes Ratio left: right')
hold on
bar(entry,'b');

hold off  

pro2 = Protein(13);
figure
grid on
axis([-1 149 -10 20])
T = strcat('Slope ratio of Protein- ', pro2);
title(T);
xlabel('Residues')
ylabel('Slopes Ratio left : right')
hold on
bar(secpro,'r');
hold off  





