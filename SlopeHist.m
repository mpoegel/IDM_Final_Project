%% Plot a histogram of the slopes

function [ h ] = SlopeHist( L, R, Ctitle )
% Adapted from Dr. Bennett's ClassHist function
% Inputs 
% L = left slopes
% R = right slopes 
% Ctitle is title of the plot
% Output h is the handle to the figure

h=figure;
hold on;

[n1, xout1] = hist(L,30);
[n2, xout2] = hist(R,xout1);
bar(xout2,[n1',n2']);
title(Ctitle);
xlabel('Slope');
ylabel('Count');
legend('Left','Right');
hold off


end

