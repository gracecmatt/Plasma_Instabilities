%%% Plot output variables from Mio

%May need to switch sign of parameter weights to put on the same scale
w= -w;

%% Plot the eigenvalues of C on a log plot
fig = figure;
semilogy(1:Nparams,evalues,'.-b','MarkerSize',30)
title(['Eigenvalues of C (N = ' int2str(N) ')'],'Interpreter','latex','Fontsize',16,'FontWeight','bold','Position',[12.5 180 0])
xlim([0,Nparams+1])
xticks(1:Nparams)
set(get(gca,'Title'),'Units','Normalized','Position',[.45,1.04])
set(fig,'PaperUnits','inches','PaperSize',[10 8])
%hgexport(fig, ['Evalues_Dispersion_' num2str(xu(2)) '_' int2str(N) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
%hgexport(fig, ['Evalues_Dispersion_' num2str(xu(2)) '_' int2str(N) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');
%hgexport(fig, ['Evalues 5% change' int2str(Nparams) '_' int2str(N) '.fig'], hgexport('factorystyle'), 'Format', 'fig');
 
%% Plot the weight vector
fig2 = figure;
plot(1:Nparams,w,'.-b','MarkerSize',30)
xlim([0,Nparams+1])
ylim([-1,1])
title(['Weight Vector (N = ' int2str(N) ')'],'Interpreter','latex','Fontsize',16,'FontWeight','bold','Position',[12.5 1.05 0])
xlabel('Parameters','Interpreter','latex','Fontsize',14)
ylabel('Parameter Weights','Interpreter','latex','Fontsize',14)
set(get(gca,'Title'),'Units','Normalized','Position',[.45,1.04])
set(gca,'TickLabelInterpreter','latex');
xticks(1:Nparams)
% Include Latex Parameter Ticks - uses alpha, k, mu, sigma
%(1) p=5
%set(gca,'XTickLabel',{'$\alpha$','$k$','$\mu$','$\sigma^2$','$\beta$'})
%(2) p=4 with beta
%set(gca,'XTickLabel',{'$\alpha$','$k$','$\sigma^2$','$\beta$'})
%(3) p=4 with mu
%set(gca,'XTickLabel',{'$\beta_1$' '$\sigma^2_1$' '$\mu_1$' '$\beta_2$' '$\sigma^2_2$' '$\mu_2$' '$k$'})
%(4) p=5
%set(gca,'XTickLabel',{'$k$','$\mu$','$\sigma^2$'})
%set(gca,'XTickLabel',{'$\alpha$','$\mu$','$\sigma^2$'})
%(5) p=6 for BiMax
set(gca,'XTickLabel',{'$k$' '$\theta$' '$\mu$'})

set(fig2,'PaperUnits','inches','PaperSize',[10 8])
%hgexport(fig2, ['WV_Dispersion_' num2str(xu(2)) '_' int2str(N) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
%hgexport(fig2, ['WV_Dispersion_' num2str(xu(2))  '_' int2str(N) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');
%hgexport(fig2, ['WV 5%change_' int2str(Nparams) '_' int2str(N) '.fig'], hgexport('factorystyle'), 'Format', 'fig');
    
%% Sufficient summary plot
fig3 = figure;
plot(Xs*w,growth,'ko');
title(['Sufficient Summary Plot (N = ' int2str(N) ')'],'Interpreter','latex','Fontsize',16,'FontWeight','bold','Position',[-.15 1130 0])
xlabel('$w^T p_j$','Interpreter','latex','FontSize',14)
ylabel('Growth Rate','Interpreter','latex','FontSize',14)
set(get(gca,'Title'),'Units','Normalized','Position',[.45,1.04])
%xlim([-2,2])
maxgrowth = max(growth); mingrowth = min(growth);
%ylim([0.9*mingrowth, 1.1*maxgrowth])
%ylim([1.1*mingrowth, 0.9*maxgrowth])
%ylim([1.1*mingrowth, 0.25])
axis square;
grid on;
set(fig3,'PaperUnits','inches','PaperSize',[10 8])
%hgexport(fig3, ['SSP_Dipsersion_' num2str(xu(2)) '_' int2str(N) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
%hgexport(fig3, ['SSP_Dispersion_' num2str(xu(2)) '_' int2str(N) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');

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
fig4 = figure;
plot(A(:,1), A(:,2), 'b');
hold on
plot(Xs*w,growth,'ko');
title(['Sufficient Summary Plot (N = ' int2str(N) ')'],'Interpreter','latex','Fontsize',16,'FontWeight','bold','Position',[-.15 1130 0])
xlabel('$w^T p_j$','Interpreter','latex','FontSize',14)
ylabel('Growth Rate','Interpreter','latex','FontSize',14)
set(get(gca,'Title'),'Units','Normalized','Position',[.45,1.04])
xlim([-2,2])
maxgrowth = max(growth); mingrowth = min(growth);
%ylim([0.9*mingrowth, 1.1*maxgrowth])
%ylim([1.1*mingrowth, 0.9*maxgrowth])
%ylim([1.1*mingrowth, 0.25])
axis square;
grid on;
set(fig4,'PaperUnits','inches','PaperSize',[10 8])
%hgexport(fig4, ['SSP_Polyfit_Dispersion_' num2str(xu(2)) '_' int2str(N) '_' int2str(deg) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
%hgexport(fig4, ['SSP_Polyfit_Dispersion_' num2str(xu(2)) '_' int2str(N) '_' int2str(deg) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');

%% Plot polynomial approximation Error
fig5 = figure;
plot(Xs*w,error,'ob')
title(['Polynomial Appx Error (degree = ',num2str(deg), ')'],'Interpreter','latex','FontSize',14)
xlabel('$w^T p_j$','Interpreter','latex','FontSize',14)
ylabel('Error','Interpreter','latex','FontSize',14)
set(get(gca,'Title'),'Units','Normalized','Position',[.45,1.04])
axis square;
grid on;
%set(fig5,'PaperUnits','inches','PaperSize',[10 8])
%hgexport(fig5, ['PolyError_Dispersion_BiMax_' int2str(100*var) '_' int2str(N) '_' int2str(deg) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
%hgexport(fig5, ['PolyError_Dispersion_BiMax_' int2str(100*var) '_' int2str(N) '_' int2str(deg) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');


%% Subplots of Eigenvalues, Weight vector, and SSP with appx

fig6 = figure('Position', [100, 100, 900, 450]);
subplot(1,3,1), semilogy(1:Nparams,evalues,'.-b','MarkerSize',30)
xlim([0,Nparams+1])
%title(['Eigenvalues of C (N = ' int2str(N) ')'],'Interpreter','latex','Fontsize',16,'FontWeight','bold','Position',[12.5 180 0])
title('Eigenvalues of C','Interpreter','latex','Fontsize',16,'FontWeight','bold','Position',[12.5 180 0])
xlabel('$\ell$','Interpreter','latex','Fontsize',14)
ylabel('$\lambda_\ell$','Interpreter','latex','Fontsize',14)
set(get(gca,'Title'),'Units','Normalized','Position',[.45,1.04])

subplot(1,3,2), plot(1:Nparams,w,'.-b','MarkerSize',30)
xlim([0,Nparams+1])
ylim([-1,1])
%title(['Weight Vector (N = ' int2str(N) ')'],'Interpreter','latex','Fontsize',16,'FontWeight','bold','Position',[12.5 1.05 0])
title('Weight Vector','Interpreter','latex','Fontsize',16,'FontWeight','bold','Position',[12.5 1.05 0])
xlabel('Parameters','Interpreter','latex','Fontsize',14)
ylabel('Parameter Weights','Interpreter','latex','Fontsize',14)
set(get(gca,'Title'),'Units','Normalized','Position',[.45,1.04])
set(gca,'TickLabelInterpreter','latex');
xticks(1:Nparams)
%set(gca,'XTickLabel',{'$\alpha$','$\mu$','$\sigma^2$'})
%set(gca,'XTickLabel',{'$k$','$\mu$','$\sigma^2$'})
%set(gca,'XTickLabel',{'$\alpha$','$k$','$\mu$','$\sigma^2$','$\beta$'})
%set(gca,'XTickLabel',{'$\beta_1$' '$\sigma^2_1$' '$\mu_1$' '$\beta_2$' '$\sigma^2_2$' '$\mu_2$' '$k$'})
set(gca,'XTickLabel',{'$k$' '$\theta$' '$\mu$'})

subplot(1,3,3), plot(A(:,1), A(:,2), 'r');
hold on
plot(Xs*w,growth,'ko');
%title(['Sufficient Summary Plot (N = ' int2str(N) ')'],'Interpreter','latex','Fontsize',16,'FontWeight','bold','Position',[-.15 1130 0])
title('Sufficient Summary Plot','Interpreter','latex','Fontsize',16,'FontWeight','bold','Position',[-.15 1130 0])
xlabel('$w^T p_j$','Interpreter','latex','FontSize',14)
ylabel('Growth Rate','Interpreter','latex','FontSize',14)
set(get(gca,'Title'),'Units','Normalized','Position',[.45,1.04])
%xlim([-2,2])
maxgrowth = max(growth); mingrowth = min(growth);
%ylim([0.9*mingrowth, 1.1*maxgrowth])
%ylim([1.1*mingrowth, 1.1*maxgrowth])
%ylim([1.1*mingrowth, 0.25])
%axis square;
grid on;

set(fig6,'PaperUnits','inches','PaperSize',[11 8])
hgexport(fig6, ['EigWVSSPfit_Dispersion_Kappa_' int2str(100*var) '_' int2str(N) '_' int2str(deg) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
%hgexport(fig6, ['EigWVSSPfit_Dispersion_BiMax_global_' int2str(N) '_' int2str(deg) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
hgexport(fig6, ['EigWVSSPfit_Dispersion_Kappa_' int2str(100*var) '_' int2str(N) '_' int2str(deg) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');
%hgexport(fig6, ['EigWVSSPfit_Dispersion_BiMax_global_' int2str(N) '_' int2str(deg) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');

%hgexport(fig6, ['EigWVSSPfit_Dispersion_BiMax_' int2str(100*var) '_' int2str(N) '_' int2str(deg) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
%hgexport(fig6, ['EigWVSSPfit_Dispersion_BiMax_' int2str(100*var) '_' int2str(N) '_' int2str(deg) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');
