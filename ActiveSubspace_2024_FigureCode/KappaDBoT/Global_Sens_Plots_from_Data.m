%% Plot output variables
% updated positions, fonts, sizing, and removed unused plots

%May need to switch sign of parameter weights to put on the same scale
w= -w;

%Compute and plot polynomial approximation & errors
deg = 2; %Order of polynomial approximation
p = polyfit(Xs*w,growth,deg);
error = abs(growth - polyval(p,Xs*w));
l2err = error'*error

%Plot polynomial approximation
%Increase grid width by a factor of gw
%NO CHANGE: gw = 1;
gw = 1.25;
minx = min(Xs*w);
maxx = max(Xs*w);
gridx = [gw*minx; Xs*w; gw*maxx];
polygrid = polyval(p,gridx);
A = [gridx, polygrid];
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
% kappa = 2 color: '#dc143c'
% kappa = 4 color: '#2f5ada'
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
legend('Quadratic Fit','Data','FontSize',12,'Location','NorthWest','Box','off');
maxgrowth = max(growth); mingrowth = min(growth);
grid on;

txt = ['\textbf{Kappa Bump-on-Tail ($\kappa\!=\!',int2str(kappa),'$), ',int2str(var*100),'\% Variation}'];
sgtitle(txt,'FontWeight','bold','Interpreter','latex','FontSize',22);

%%
set(gcf, 'PaperPosition', [0 0 9 6]); %Position the plot further to the upper-left conder
set(gcf, 'PaperSize', [9 6]); % Extends the plot to fill the entire paper
saveas(gcf, ['Figs/EigWVSSPfit_Dispersion_BiKap' int2str(kappa) '_' int2str(Nparams) '_' int2str(N) '_' int2str(100*var) '_' int2str(deg) '.pdf'])
saveas(gcf, ['Figs/EigWVSSPfit_Dispersion_BiKap' int2str(kappa) '_' int2str(Nparams) '_' int2str(N) '_' int2str(100*var) '_' int2str(deg) '.fig'])