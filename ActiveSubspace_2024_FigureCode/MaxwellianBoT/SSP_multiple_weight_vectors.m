%% Plotting first two weight vectors
%clear all; close all; 

% Load workspace 
%load('Dispersion_Rate_BiMax_P7_N512_0.05data_par.mat')
%load('Dispersion_Rate_V2_P3_N512_0.25data_par2.mat')
deg = 3; %???
w_1_1var = w; 
w_2_1var = w2; 
growth_1var = growth; 
Xs_25var = Xs; 
eta(1) = (evalues(1))/sum(evalues);
eta(2) = (evalues(1)+evalues(2))/sum(evalues);


set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure; 
scatter(Xs*w_1_1var, Xs*w_2_1var, 60, growth_1var, 'o','filled');
colormap(jet);
c1 = colorbar;
xlabel('$w_1^Tp_j$','Fontsize',14)
ylabel('$w_2^Tp_j$','Fontsize',14)
ylabel(c1,'Growth Rate','FontSize',14,'Rotation',270,'Interpreter','latex');
c1.Label.Position(1) = 4;
%title('15% variation')
% title('$25\%$ variation','Interpreter','latex','Fontsize',16,'FontWeight','bold')
%title('Global parameter variation','Interpreter','latex','Fontsize',16,'FontWeight','bold')
txt = ['Maxwellian BoT, ',int2str(var*100),'\% Variation'];
title(txt,'FontSize',16);

set(gcf, 'PaperPosition', [0 0 6 4.5]); %Position the plot further to the upper-left conder
set(gcf, 'PaperSize', [6 4.5]); % Extends the plot to fill the entire paper
saveas(gcf, ['Figs/Eig2fit_Dispersion_BiMax_' int2str(Nparams) '_' int2str(N) '_' int2str(100*var) '_' int2str(deg) '.pdf'])
saveas(gcf, ['Figs/Eig2fit_Dispersion_BiMax_' int2str(Nparams) '_' int2str(N) '_' int2str(100*var) '_' int2str(deg) '.fig'])


%% 3D plot of growth rate vs first two weight vectors
% figure;
% tri = delaunay(Xs*w_1_1var, Xs*w_2_1var);
% m = trimesh(tri, Xs*w_1_1var, Xs*w_2_1var, growth_1var);
% axis vis3d
% shading interp
% colorbar EastOutside
% xlabel('$w_1^Tp_j$', 'Interpreter','latex','Fontsize',16,'FontWeight','bold')
% ylabel('$w_2^Tp_j$', 'Interpreter','latex','Fontsize',16,'FontWeight','bold')
% %ylabel(c1,'Growth Rate','FontSize',16,'Rotation',270);
% %c1.Label.Position(1) = 5;
% 
% figure;
% s = trisurf(tri, Xs*w_1_1var, Xs*w_2_1var, growth_1var);
% axis vis3d
% shading interp
% colorbar EastOutside
% xlabel('$w_1^Tp_j$', 'Interpreter','latex','Fontsize',16,'FontWeight','bold')
% ylabel('$w_2^Tp_j$', 'Interpreter','latex','Fontsize',16,'FontWeight','bold')
% %ylabel(c1,'Growth Rate','FontSize',16,'Rotation',270);
% %c1.Label.Position(1) = 5;
% %colormap default;
% %cb = colorbar;
% %ylabel(cb,'$\gamma(k, \sigma^2)$','Interpreter','latex','FontSize',14,'Rotation',90)

x = Xs*w_1_1var;
y = Xs*w_2_1var;
z = growth_1var;

M = 1000;
xi = linspace(min(x),max(x),M) ;
yi = linspace(min(y),max(y),M) ;
[Xi,Yi] = meshgrid(xi,yi) ;
% Zi = griddata(x,y,z,Xi,Yi) ;
% figure;
% mesh(Xi,Yi,Zi);
% %surf(Xi,Yi,Zi);
% shading interp
% colorbar EastOutside
% xlabel('$w_1^Tp_j$', 'Interpreter','latex','Fontsize',16,'FontWeight','bold')
% ylabel('$w_2^Tp_j$', 'Interpreter','latex','Fontsize',16,'FontWeight','bold')

%% Nonlinear curve fit for 2D active subspace 
x = Xs*w_1_1var;
y = Xs*w_2_1var;
z = growth_1var;


f12 = fit( [x, y], z, 'poly12');

%Plot polynomial approximation
fig = figure;
plot(f12, [x,y], z);
colormap(jet);
c1 = colorbar;
title(['2D Sufficient Summary Plot (N = ' int2str(N) ')'],'Interpreter','latex','Fontsize',16,'FontWeight','bold')
xlabel('$w_1^T p_j$','Interpreter','latex','FontSize',16)
ylabel('$w_2^T p_j$','Interpreter','latex','FontSize',16)
ylabel(c1,'Growth Rate','FontSize',16,'Rotation',270);
c1.Label.Position(1) = 4;
set(get(gca,'Title'),'Units','Normalized')


%maxgrowth = max(z); mingrowth = min(z);
%zlim([0.8*maxgrowth, 1.2*maxgrowth])
%axis square;
grid on;
set(fig,'PaperUnits','inches','PaperSize',[11 8])
%hgexport(fig, ['SSP2_Dipsersion_' int2str(N) '_global.eps'], hgexport('factorystyle'), 'Format', 'eps');
%hgexport(fig, ['SSP2_Dipsersion_' int2str(N) '_global.pdf'], hgexport('factorystyle'), 'Format', 'pdf');

%Compute and plot polynomial approximation & errors
error = abs(growth - f12(x,y));
l2err = error'*error


%% Plotting first three weight vectors
% figure; 
% 
% scatter3(Xs*w_1_1var, Xs*w_2_1var, Xs*w_3_1var, 60, growth_1var, 'o', 'filled');
% colormap(jet);
% c2 = colorbar;
% xlabel('$w_1^Tp_j$', 'Interpreter','latex','Fontsize',16,'FontWeight','bold')
% ylabel('$w_2^Tp_j$', 'Interpreter','latex','Fontsize',16,'FontWeight','bold')
% ylabel(c2,'Growth Rate','FontSize',16,'Rotation',270);
% zlabel('$w_3^Tp_j$', 'Interpreter','latex','Fontsize',16,'FontWeight','bold')
% c2.Label.Position(1) = 5;
% title('V^2-Two Stream: 25% variation')
% ax = gca; 
% ax.FontSize = 16; 
% hold on
% 
% 
% 
% %set(gcf,'Units','Normalized','OuterPosition',[0 0 1 1])
% saveas(gcf, 'SSP_V2_25%var_3D', 'epsc')

