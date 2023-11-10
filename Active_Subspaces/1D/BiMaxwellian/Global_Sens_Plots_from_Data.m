%% Plot output variables from Mio
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

%% Subplots of Eigenvalues, Weight vector, and SSP with appx
close all

set(0,'DefaultAxesFontName','Arial')
fig6 = figure('Position', [100, 100, 900, 450]);

% SUBPLOT 1: EIGENVALUES
pos = [0.1 0.1 0.225 0.74];
subplot('Position',pos), semilogy(1:Nparams,evalues,'.-b','MarkerSize',30)
ax = gca; ax.FontSize = 10; %use for setting tick font size
xlim([0,Nparams+1])
title('Eigenvalues of C','Fontsize',16);
xlabel('$\ell$','Interpreter','latex','Fontsize',15)
ylabel('$\lambda_\ell$','Interpreter','latex','Fontsize',15)

% SUBPLOT 2: WEIGHT VECTOR
pos = [0.4 0.1 0.225 0.74];
subplot('Position',pos), plot(1:Nparams,w,'.-b','MarkerSize',30)
ax = gca; ax.FontSize = 10; %use for setting tick font size
    xaxisproperties=get(gca, 'XAxis');
    xaxisproperties.TickLabelInterpreter = 'latex'; %latex for x-axis only
    xaxisproperties.FontSize = 12; %enlarge latex font
xlim([0,Nparams+1]); % ylim([-1,1])
title('Weight Vector','Fontsize',16);
xlabel('Parameters','Fontsize',14)
ylabel('Parameter Weights','Fontsize',14)
xticks(1:Nparams)
set(gca,'XTickLabel',{'$k$' '$\sigma_1$' '$\sigma_2$' '$\mu_1$' '$\mu_2$' '$\beta$'})

% SUBPLOT 3: SUFFICIENT SUMMARY
pos = [0.7 0.1 0.225 0.74];
subplot('Position',pos);
plot(A(:,1), A(:,2), 'r'); hold on
plot(Xs*w,growth,'ko');
ax = gca; ax.FontSize = 10; %use for setting tick font size
title('Sufficient Summary Plot','Fontsize',16);
xlabel('$w^T p_j$','Interpreter','latex','FontSize',15);
xlabel('$\gamma(p)$','Interpreter','latex','FontSize',15);
legend('Quadratic Fit','Data','FontSize',9.5,'Location','NorthWest','Box','off');
maxgrowth = max(growth); mingrowth = min(growth);
grid on;

set(fig6,'PaperUnits','inches','PaperSize',[11 8.5])
txt = ['BiMaxwellian, ',int2str(var*100),'% Variation'];
sgtitle(txt,'Fontsize',23,'FontWeight','bold');

hgexport(fig6, ['Figs\EigWVSSPfit_Dispersion_BiMax_' int2str(100*var) '_' int2str(N) '_' int2str(deg) '.eps']);
hgexport(fig6, ['Figs\EigWVSSPfit_Dispersion_BiMax_' int2str(100*var) '_' int2str(N) '_' int2str(deg) '.pdf']);