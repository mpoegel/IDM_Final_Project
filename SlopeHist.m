
function [ h ] = SlopeHist( L, R, Ctitle )
%% Plot class histograms and threshold line
% Inputs 
% Cp = Scalar projections of Class 1
% Cm = Scalar projections of Class -1 
% t= threhold
% err=error in percentage
% Title is title of the plot
% Output h is the handle to the figure
%% Start the Figure
h=figure;
hold on;

%%Calculate the binsizes
min_val = min([L;R]);
max_val = max([L;R]);
binsize=(max_val-min_val)/20;

%% Plot the histogram  

[n1, xout1] = hist(L);
[n2, xout2] = hist(R,xout1);
bar(xout2,[n1',n2']);
title(Ctitle);
xlabel('scalar projection');
ylabel('count');
legend('Left','Right');
hold off


end

