%% Plot output variables
%May need to switch sign of parameter weights to put on the same scale
w= -w;

%Compute quadratic polynomial fit 
[fitted_curve1,gof1] = fit(Xs*w,growth,'poly2');
% Save the coeffiecient values for p1, p2, and p3 in a vector
p1 = coeffvalues(fitted_curve1);
poly2error = [gof1.sse; gof1.rmse; gof1.rsquare; gof1.adjrsquare];

%Compute atan fit
x0 = -[1 1 1 1]; 
fitfun = fittype( @(a,b,c,d,x) d*atan(b*(x+a))+c );
[fitted_curve2,gof2] = fit(Xs*w,growth,fitfun,'StartPoint',x0);
% Save the coeffiecient values for a,b,c and d in a vector
p2 = coeffvalues(fitted_curve2);
atanerror = [gof2.sse; gof2.rmse; gof2.rsquare; gof2.adjrsquare];

error = table(poly2error,atanerror,'VariableNames',["Poly2 Error","Arctan Error"],'RowNames',["sse","rmse","rsquare","adjrsquare"])

% change type of error to compare by picking element 1-4
% if error.("Poly2 Error")(1)<error.("Arctan Error")(1) %choose quadratic fit
%     p = p1;
%     fitcurve = fitted_curve1;
%     deg = 2;
% elseif error.("Poly2 Error")(1)>=error.("Arctan Error")(1) %choose arctan fit
    p = p2;
    fitcurve = fitted_curve2;
    deg = 4;
% end
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
if kappa ==2; color = '#dc143c'; else; color = '#2f5ada'; end
% colorSet = {'#dc143c','#cb2888','#9932cc','#6c4ce2','#2f5ada'};%red->blue
% kappa = 2 color: '#dc143c'; kappa = 4 color: '#2f5ada'
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% SUBPLOT 1: EIGENVALUES
pos = [0.1 0.1 0.225 0.74];
subplot('Position',pos), semilogy(1:Nparams,evalues,'.-','Color',color,'MarkerSize',30,'linewidth',2)
ax = gca; ax.FontSize = 10; %use for setting tick font size
xlim([0,Nparams+1]);
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
xlim([0,Nparams+1]);  ylim([-1,1]);
set(gca, 'YGrid', 'on', 'XGrid', 'off', 'YMinorGrid', 'off');
title('Weight Vector','Fontsize',16);
xlabel('Parameters','Fontsize',14)
ylabel('Parameter Weights','Fontsize',14)
xticks(1:Nparams)
set(gca,'XTickLabel',{'$k$' '$\sigma_1$' '$\sigma_2$' '$\mu$' '$v_0$' '$\beta$' '$\kappa$'})

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

txt = ['\textbf{Kappa Bump-on-Tail ($\kappa\!=\!',int2str(kappa),'$), ',int2str(var*100),'\% Variation}'];
sgtitle(txt,'FontWeight','bold','Interpreter','latex','FontSize',22);

set(gcf, 'PaperPosition', [0 0 9 6]); %Position the plot further to the upper-left conder
set(gcf, 'PaperSize', [9 6]); % Extends the plot to fill the entire paper
saveas(gcf, ['Figs/EigWVSSPfit_Dispersion_KapBoT' int2str(kappa) '_' int2str(Nparams) '_' int2str(N) '_' int2str(100*var) '_' int2str(deg) '.pdf'])
saveas(gcf, ['Figs/EigWVSSPfit_Dispersion_KapBoT' int2str(kappa) '_' int2str(Nparams) '_' int2str(N) '_' int2str(100*var) '_' int2str(deg) '.svg'])
saveas(gcf, ['Figs/EigWVSSPfit_Dispersion_KapBoT' int2str(kappa) '_' int2str(Nparams) '_' int2str(N) '_' int2str(100*var) '_' int2str(deg) '.fig'])

save(['Data/Fits_KappaBoT' int2str(kappa) '_var' int2str(var*100) '_data.mat'],"y","fitcurve")
