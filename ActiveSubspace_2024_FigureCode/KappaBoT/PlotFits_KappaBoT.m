clear; clc;

data = load('C:\Users\15134\Documents\Research\Local_Github\ActiveSubspace_2024_FigureCode\KappaBoT\Data\Fits_KappaBoT2_var1_data.mat');
y1 = data.y;
fitcurve1 = data.fitcurve;
A1 = [y1,fitcurve1(y1)];
[~, order] = sort(A1(:,1));
A1 = A1(order,:);

data = load('C:\Users\15134\Documents\Research\Local_Github\ActiveSubspace_2024_FigureCode\KappaBoT\Data\Fits_KappaBoT2_var5_data.mat');
y5 = data.y;
fitcurve5 = data.fitcurve;
A5 = [y5,fitcurve5(y5)];
[~, order] = sort(A5(:,1));
A5 = A5(order,:);

data = load('C:\Users\15134\Documents\Research\Local_Github\ActiveSubspace_2024_FigureCode\KappaBoT\Data\Fits_KappaBoT2_var10_data.mat');
y10 = data.y;
fitcurve10 = data.fitcurve;
A10 = [y10,fitcurve10(y10)];
[~, order] = sort(A10(:,1));
A10 = A10(order,:);

data = load('C:\Users\15134\Documents\Research\Local_Github\ActiveSubspace_2024_FigureCode\KappaBoT\Data\Fits_KappaBoT2_var15_data.mat');
y15 = data.y;
fitcurve15 = data.fitcurve;
A15 = [y15,fitcurve15(y15)];
[~, order] = sort(A15(:,1));
A15 = A15(order,:);

data = load('C:\Users\15134\Documents\Research\Local_Github\ActiveSubspace_2024_FigureCode\KappaBoT\Data\Fits_KappaBoT2_var20_data.mat');
y20 = data.y;
fitcurve20 = data.fitcurve;
A20 = [y20,fitcurve20(y20)];
[~, order] = sort(A20(:,1));
A20 = A20(order,:);

data = load('C:\Users\15134\Documents\Research\Local_Github\ActiveSubspace_2024_FigureCode\KappaBoT\Data\Fits_KappaBoT2_var25_data.mat');
y25 = data.y;
fitcurve25 = data.fitcurve;
A25 = [y25,fitcurve25(y25)];
[~, order] = sort(A25(:,1));
A25 = A25(order,:);

% data = load('C:\Users\15134\Documents\Research\Local_Github\ActiveSubspace_2024_FigureCode\KappaBoT\Data\Fits_KappaBoT2_var50_data.mat');
% y50 = data.y;
% fitcurve50 = data.fitcurve;
% A50 = [y50,fitcurve50(y50)];
% [~, order] = sort(A50(:,1));
% A50 = A50(order,:);

%% Plot
close all;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure
plot(A1(:,1),A1(:,2),'linewidth',2); hold on
plot(A5(:,1),A5(:,2),'linewidth',2); 
plot(A10(:,1),A10(:,2),'linewidth',2); 
plot(A15(:,1),A15(:,2),'linewidth',2); 
plot(A20(:,1),A20(:,2),'linewidth',2); 
plot(A25(:,1),A25(:,2),'linewidth',2); 
% plot(A50(:,1),A50(:,2),'linewidth',2); 
ylabel("$\gamma(y)$ fit function",'Fontsize',14);
xlabel("$y=\omega^T p$",'Fontsize',14); xlim([-1.75,1.75])
title("Function Fits With Increasing Variation",'Fontsize',16)
legend("1\% Variation","5\% Variation","10\% Variation","15\% Variation","20\% Variation","25\% Variation",'Location','Best','Fontsize',12)
grid on