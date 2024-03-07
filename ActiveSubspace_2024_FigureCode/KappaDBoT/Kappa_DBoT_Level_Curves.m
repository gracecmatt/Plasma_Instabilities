clear; clc;
kappa = 2;
v0_grid = linspace(2,6,60);
beta_grid = linspace(0.5,1,60);
parpool(10);

v0_length = length(v0_grid); beta_length = length(beta_grid);
gamma_grid = zeros(v0_length,beta_length);

tic;
parfor ii=1:v0_length
    for jj=1:beta_length

        v0 = v0_grid(ii);
        beta = beta_grid(jj);

        params = [0.5; 1; 0.7; 1; v0; beta; kappa];

        init_guess = Vlasov_1D_linearized_Steve_v4_Kappa(params(1),params(2),params(3),0,params(5),params(6),params(7));
        xi_guess = init_guess/(params(2)*params(1)); % shifted & scaled

        omega = BiKappa_Disp_Using_Xie(params(1)*params(2),1,params(3)/params(2),0,params(5)/params(2),params(6),params(7),xi_guess)*params(2)*params(1);% + params(4)*params(1);
        gamma_grid(ii,jj) = imag(omega);

    end
end
delete(gcp('nocreate')); toc

%%

figure;
surf(v0_grid, beta_grid, gamma_grid);
view(0, 90);
cb = colorbar;
colormap(bluewhitered);
xlabel('$v_0$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$\beta$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel(cb, 'Growth Rate $\gamma$', 'FontSize', 16, 'Interpreter', 'latex');
ax = gca; ax.FontSize = 20;
xlim([min(v0_grid), max(v0_grid)]);
ylim([min(beta_grid), max(beta_grid)]);
title('Kappa DBoT Level Curves: $v_0$, $\beta$ vs. $\gamma$', ...
    'Interpreter', 'latex', 'FontSize', 20);
params = [0.5; 1; 0.7; 1; -999; -999; kappa];
subtitle(['$k=',num2str(params(1)),',\sigma_1=',num2str(params(2)), ...
    ',\sigma_2=',num2str(params(3)),',\mu=',num2str(params(4)),',\kappa=',num2str(params(7)),'$']);

set(gcf, 'PaperPosition', [0 0 9 6]); %Position the plot further to the upper-left conder
set(gcf, 'PaperSize', [9 6]); % Extends the plot to fill the entire paper
saveas(gcf, 'Figs/Kappa_DBoT_Level_Curves_v0_beta_vs_gamma.pdf')
saveas(gcf, 'Figs/Kappa_DBoT_Level_Curves_v0_beta_vs_gamma.fig')

save('Data/Kappa_DBoT_Level_Curves_v0_beta_vs_gamma')
