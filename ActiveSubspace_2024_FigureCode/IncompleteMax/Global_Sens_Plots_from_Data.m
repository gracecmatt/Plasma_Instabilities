%% Plot output variables
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
% close all
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
set(gca,'XTickLabel',{'$k$' '$\sigma$' '$\mu$' '$\nu$'})
 
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
% ylim([0.9*mingrowth, 1.1*maxgrowth])
% ylim([-0.35,0.05])
grid on;

txt = ['\textbf{Incomplete Maxwellian, ',int2str(var*100),'\% Variation}'];
sgtitle(txt,'FontWeight','bold','FontSize',22);

set(gcf, 'PaperPosition', [0 0 9 6]); %Position the plot further to the upper-left conder
set(gcf, 'PaperSize', [9 6]); % Extends the plot to fill the entire paper
saveas(gcf, ['Figs/EigWVSSPfit_Dispersion_IMax_' int2str(Nparams) '_' int2str(N) '_' int2str(100*var) '_' int2str(deg) '.pdf'])
saveas(gcf, ['Figs/EigWVSSPfit_Dispersion_IMax_' int2str(Nparams) '_' int2str(N) '_' int2str(100*var) '_' int2str(deg) '.fig'])