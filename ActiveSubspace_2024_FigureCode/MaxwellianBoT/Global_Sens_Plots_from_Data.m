%% Plot output variables
%May need to switch sign of parameter weights to put on the same scale
w= -w;

% % %Compute quadratic polynomial fit 
% [fitted_curve1,gof1] = fit(Xs*w,growth,'poly3');
% % Save the coeffiecient values for p1, p2, and p3 in a vector
% p1 = coeffvalues(fitted_curve1);
% poly2error = [gof1.sse; gof1.rmse; gof1.rsquare; gof1.adjrsquare];

%Compute atan fit with +x0 
fitfun = fittype( @(a,b,c,d,x) a*atan(b*x+c)+d );
x0_1 = [1 1 1 1];
[fitted_curve1,gof] = fit(Xs*w,growth,fitfun,'StartPoint',x0_1);
% Save the coeffiecient values for a,b,c and d in a vector
p1 = coeffvalues(fitted_curve1);
atanerror1 = [gof.sse; gof.rmse; gof.rsquare; gof.adjrsquare];

%Compute atan fit with -x0
fitfun = fittype( @(a,b,c,d,x) a*atan(b*x+c)+d );
x0_2 = [-1 1 1 1]; 
[fitted_curve2,gof2] = fit(Xs*w,growth,fitfun,'StartPoint',x0_2);
% Save the coeffiecient values for a,b,c and d in a vector
p2 = coeffvalues(fitted_curve2);
atanerror2 = [gof2.sse; gof2.rmse; gof2.rsquare; gof2.adjrsquare];

error = table(atanerror1,atanerror2,'VariableNames',["Arctan Error +","Arctan Error -"],'RowNames',["sse","rmse","rsquare","adjrsquare"])

% change type of error to compare by picking element 1-4
if error.("Arctan Error +")(1)<error.("Arctan Error -")(1) %choose arctan fit +
    p = p1;
    fitcurve = fitted_curve1;
    deg = 4;
elseif error.("Arctan Error +")(1)>=error.("Arctan Error -")(1) %choose arctan fit -
    p = p2;
    fitcurve = fitted_curve2;
    deg = 4;
end
gw = 1.25;
y = Xs*w;
minx = min(Xs*w);
maxx = max(Xs*w);
gridx = [gw*minx; Xs*w; gw*maxx];
A = [gridx, fitcurve(gridx)];
[temp, order] = sort(A(:,1));
A = A(order,:);

for i=1:length(evalues)
    if evalues(i)<10^(-16); evalues(i)=0; end
end

%% Subplots of Eigenvalues, Weight vector, and SSP with appx
close all
fig6 = figure('Position', [100, 100, 900, 450]);
% colorSet = {'#dc143c','#cb2888','#9932cc','#6c4ce2','#2f5ada'};%red->blue
color = '#2f5ada';
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% SUBPLOT 1: EIGENVALUES
pos = [0.1 0.1 0.225 0.74];
subplot('Position',pos), semilogy(1:Nparams,evalues,'.-','Color',color,'MarkerSize',30,'linewidth',2)
ax = gca; ax.FontSize = 10; %use for setting tick font size
xlim([0,Nparams+1])
set(gca, 'YGrid', 'on', 'XGrid', 'off', 'YMinorGrid', 'off');
title('Eigenvalues of C','Fontsize',16);
xlabel('$\ell$','Fontsize',14)
ylabel('$\lambda_\ell$','Fontsize',14)

% SUBPLOT 2: WEIGHT VECTOR
pos = [0.4 0.1 0.225 0.74];
subplot('Position',pos), plot(1:Nparams,w,'.-','Color',color,'MarkerSize',30,'linewidth',2)
ax = gca; ax.FontSize = 10; %use for setting tick font size
    xaxisproperties=get(gca, 'XAxis');
    xaxisproperties.TickLabelInterpreter = 'latex'; %latex for x-axis only
    xaxisproperties.FontSize = 12; %enlarge latex font
xlim([0,Nparams+1]);  ylim([-1,1])
set(gca, 'YGrid', 'on', 'XGrid', 'off', 'YMinorGrid', 'off');
title('Weight Vector','Fontsize',16);
xlabel('Parameters','Fontsize',14)
ylabel('Parameter Weights','Fontsize',14)
xticks(1:Nparams)
% set(gca,'XTickLabel',{'$k$' '$\sigma_1$' '$\sigma_2$' '$\mu_1$' '$\mu_2$' '$\beta$'})
set(gca,'XTickLabel',{'$k$' '$\mu_1$' '$\mu_2$' '$\sigma_1$' '$\sigma_2$' '$\beta$'}) %switch the order to match paper
 
% SUBPLOT 3: SUFFICIENT SUMMARY
pos = [0.7 0.1 0.225 0.74];
subplot('Position',pos);
plot(A(:,1), A(:,2), 'Color',color,'linewidth',1.5); hold on
plot(Xs*w,growth,'ko');
ax = gca; ax.FontSize = 10; %use for setting tick font size
title('Sufficient Summary Plot','Fontsize',16);
xlabel('$w^T p_j$','FontSize',14);
ylabel('$\gamma(p)$','FontSize',14);
legend('Fit','Data','FontSize',12,'Location','NorthWest','Box','off');
maxgrowth = max(growth); mingrowth = min(growth);
grid on;

txt = ['\textbf{Maxwellian Bump-on-Tail, ',int2str(var*100),'\% Variation}'];
sgtitle(txt,'FontWeight','bold','FontSize',22);

set(gcf, 'PaperPosition', [0 0 9 6]); %Position the plot further to the upper-left conder
set(gcf, 'PaperSize', [9 6]); % Extends the plot to fill the entire paper
saveas(gcf, ['Figs/EigWVSSPfit_Dispersion_BiMax_' int2str(Nparams) '_' int2str(N) '_' int2str(100*var) '_' int2str(deg) '.pdf'])
saveas(gcf, ['Figs/EigWVSSPfit_Dispersion_BiMax_' int2str(Nparams) '_' int2str(N) '_' int2str(100*var) '_' int2str(deg) '.fig'])