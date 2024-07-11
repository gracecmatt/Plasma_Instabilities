%% New Plot
close all
clear;clc;

% colorSet = {'#dc143c','#cb2888','#9932cc','#6c4ce2','#2f5ada'};
% colors = {'#0d3dbf', '#4770cb','#6388d0','#8325b3','#ba567a','#f08742'};
% colors2 = {'#ffa600','#ff8531','#ff6361','#bc5090','#8a508f','#2c4875'};
colors = {'#A2142F','#D95319','#EDB120','#77AC30','#4DBEEE','#0072BD'};
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% 5% Variation around p=(0.5, 1, 1.5, 1, 4)
w1_kapBoT2 = -[0.594902088918320;0.116678223082414;-0.033205539689206;0;0.512124167353664;0.606606014478381];
w1_kapBoT4 = -[0.576095044576059;0.198267452491796;0.064766642863697;0;0.312737478867713;0.725048120482877];
w1_kapBoT6 = [-0.573082491876338;-0.233133028538152;-0.101365043031907;0;-0.237818435076973;-0.741400964175733];
w1_kapBoT8 = [-0.563717663909931;-0.249715057566427;-0.108562617404332;0;-0.205067036333643;-0.751987094946160];
w1_kapBoT50 = [-0.547777413114844;-0.282926126081842;-0.146701470345636;0;-0.117381721554880;-0.764573273616445]; 
w1_maxBoT = [-0.540042646359481;-0.285697252553735;-0.153678887639412;0;-0.099030828787275;-0.770264055009653];

K = 5;
w1 = zeros(6,K);

for i = 1:6
    w1(i,:) = [w1_kapBoT2(i),w1_kapBoT4(i),w1_kapBoT6(i),w1_kapBoT8(i),w1_maxBoT(i)];
end

plot([1,2,3,4,10],w1,'.-','MarkerSize',20,'linewidth',2); hold on;
ax = gca; ax.FontSize = 12; %use for setting tick font size
    xaxisproperties=get(gca, 'XAxis');
    xaxisproperties.TickLabelInterpreter = 'latex'; %latex for x-axis only
    xaxisproperties.FontSize = 14; %enlarge latex font
xlim([0,11.6]);  ylim([-1,0.25]);
set(gca, 'YGrid', 'on', 'XGrid', 'off', 'YMinorGrid', 'off');
colororder(colors);

txt = {'\textbf{Active Weight Vector of Bump-on-Tail Distributions}';'\textbf{at 5\% Variation around $\mathbf{p}=(0.5,1,1.5,1,4,0.9)$}'};
title(txt,'Fontsize',20);
xlabel('$\kappa$ Index','Fontsize',18)
ylabel('Parameter Weights','Fontsize',18)
xticks([1,2,3,4,10])
set(gca,'XTickLabel',{'$\kappa=2$' '$\kappa=4$' '$\kappa=6$' '$\kappa=8$' '$\kappa=\infty$'})
legend('$k$','$\sigma_1$' ,'$\sigma_2$' ,'$\mu$' ,'$v_0$', '$\beta$','FontSize',12,'Location','East');


set(gcf, 'PaperPosition', [0 0 9 3.5]); %Position the plot further to the upper-left conder
set(gcf, 'PaperSize', [9 3.5]); % Extends the plot to fill the entire paper
saveas(gcf, ['compareBoTWeightVec.pdf'])
saveas(gcf, ['compareBoTWeightVec.svg'])
saveas(gcf, ['compareBoTWeightVec.fig'])

%%
colors = {'#A2142F','#D95319','#EDB120','#77AC30','#4DBEEE','#0072BD'};
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

meangammasBOT = [-0.0248,-0.009972358891841,-9.2841e-04,0.005619048579919,0.0104,0.017952815531282,0.024397512553057,0.0281,0.0317];
kappasBOT = [2,3,4,5,6,10,20,50,1000];
meangammas = [-0.175387686622341,-0.228011315063962,-0.229680865194287,-0.229032625797779,-0.224328143711837,-0.213290210173026,-0.200061614593176,-0.189628900586748,-0.188475966245519];
kappas = [2,3,4,5,6,10,20,50,1000];


figure
semilogx(kappasBOT,meangammasBOT,'.-','MarkerSize',20,'linewidth',2); hold on;
semilogx(kappas,meangammas,'.-','MarkerSize',20,'linewidth',2); hold on;
ax = gca; ax.FontSize = 12; %use for setting tick font size
    xaxisproperties=get(gca, 'XAxis');
    xaxisproperties.TickLabelInterpreter = 'latex'; %latex for x-axis only
    xaxisproperties.FontSize = 14; %enlarge latex font
xlim([2,1001]); % ylim([-1,0.25]);
xticks([2,4,6,10,20,50,1000])
set(gca,'XTickLabel',{'$2$' '$4$' '$6$' '$10$' '$20$' '$50$' '$\infty$'})

set(gca, 'YGrid', 'on', 'XGrid', 'on', 'YMinorGrid', 'off', 'XMinorGrid', 'off');
colororder(colors([1,end]));
legend('Kappa BoT','Kappa','location','east','fontsize',14)
title('\textbf{Mean Growth/Decay Rate of Kappa Distributions}','fontsize',24)
xlabel('$\kappa$ Index','fontsize',17); ylabel('mean$(\gamma)$','fontsize',17)
set(gcf, 'PaperPosition', [0 0 9 3]); %Position the plot further to the upper-left conder
set(gcf, 'PaperSize', [9 3]); % Extends the plot to fill the entire paper
% saveas(gcf, ['compareMeanGamma.pdf'])
saveas(gcf, ['compareMeanGamma3.svg'])
% saveas(gcf, ['compareMeanGamma.fig'])