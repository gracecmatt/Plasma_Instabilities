%% Kappa
% save all the figures using data from a .mat file
% load the data, then run this script.

%% Plot Level Curves and Percent Errors
close all
txtbase = string(['$(k,\sigma,\mu)$=(',num2str(baseparams(1)),...
    ',',num2str(baseparams(2)),',',num2str(baseparams(3)),')']);
txtleg1 = {'Bohm-Gross','Xie','Shifted Xie','Shifted & Scaled Xie','Dielectric Func.'};
txtleg2 = {'Bohm-Gross','Xie','Shifted Xie','Shifted & Scaled Xie'};

newcolors = {'#4363d8','#e6194B','#3cb44b','#7E2F8E','#42d4f4'};

% plot the level curves
fig1 = figure; tiledlayout(2,2);
for i=1:Nparams
    nexttile
    plot(paramsarr(i,:), imag(omega.spectral(i,:)),'.-'); hold on
    plot(paramsarr(i,:), imag(omega.xie(i,:)),'.-');
    plot(paramsarr(i,:), imag(omega.xie_shift(i,:)),'.-');
    plot(paramsarr(i,:), imag(omega.xie_shiftscaled(i,:)),'.-');
    plot(paramsarr(i,:), imag(omega.dielectric(i,:)),'k');
        title('$\gamma('+txtparams(i)+')$ Level Curves','Interpreter','latex','FontSize',14);
        xlabel('$'+txtparams(i)+'$','Interpreter','latex','FontSize',12);
        ylabel('$\gamma('+txtparams(i)+')$','Interpreter','latex','FontSize',12);
        colororder(newcolors);
end
    leg = legend(txtleg1,'Orientation','Horizontal','FontSize',8);
    leg.Layout.Tile = 'north';
    sgtitle('Kappa $\gamma(p)$ around '+txtbase,'Interpreter','latex','FontSize',16);

fig2 = figure; tiledlayout(2,2);
for i=1:Nparams
    nexttile
    plot(paramsarr(i,:), real(omega.spectral(i,:)),'.-'); hold on
    plot(paramsarr(i,:), real(omega.xie(i,:)),'.-');
    plot(paramsarr(i,:), real(omega.xie_shift(i,:)),'.-');
    plot(paramsarr(i,:), real(omega.xie_shiftscaled(i,:)),'.-');
    plot(paramsarr(i,:), real(omega.dielectric(i,:)),'k');
        title('$\Omega('+txtparams(i)+')$ Level Curves','Interpreter','latex','FontSize',14);
        xlabel('$'+txtparams(i)+'$','Interpreter','latex','FontSize',12);
        ylabel('$\Omega('+txtparams(i)+')$','Interpreter','latex','FontSize',12);
        colororder(newcolors);
end
    leg = legend(txtleg1,'Orientation','Horizontal','FontSize',8);
    leg.Layout.Tile = 'north';
    sgtitle('Kappa $\Omega(p)$ around '+txtbase,'Interpreter','latex','FontSize',16);

% plot the relative error
fig3 = figure; tiledlayout(2,2);
for i=1:Nparams
    nexttile
    semilogy(paramsarr(i,:), errorIm.spectral(i,:),'.-'); hold on
    semilogy(paramsarr(i,:), errorIm.xie(i,:),'.-');
    semilogy(paramsarr(i,:), errorIm.xie_shift(i,:),'.-');
    semilogy(paramsarr(i,:), errorIm.xie_shiftscaled(i,:),'.-');
        title('Relative Error in $\gamma('+txtparams(i)+')$','Interpreter','latex','FontSize',14);
        xlabel('$'+txtparams(i)+'$','Interpreter','latex','FontSize',12);
        ylabel('Relative Error','Interpreter','latex','FontSize',12);
        ylim([0,max(1,max(errorIm.spectral(i,:)))]); grid on; colororder(newcolors);
end
    leg = legend(txtleg2,'Orientation','Horizontal','FontSize',9);
    leg.Layout.Tile = 'north';
    sgtitle('Kappa Error in $\gamma(p)$ around '+txtbase,'Interpreter','latex','FontSize',16);

fig4 = figure; tiledlayout(2,2);
for i=1:Nparams
    nexttile
    semilogy(paramsarr(i,:), errorRe.spectral(i,:),'.-'); hold on
    semilogy(paramsarr(i,:), errorRe.xie(i,:),'.-');
    semilogy(paramsarr(i,:), errorRe.xie_shift(i,:),'.-');
    semilogy(paramsarr(i,:), errorRe.xie_shiftscaled(i,:),'.-');
        title('Relative Error in $\Omega('+txtparams(i)+')$','Interpreter','latex','FontSize',14);
        xlabel('$'+txtparams(i)+'$','Interpreter','latex','FontSize',12);
        ylabel('Relative Error','Interpreter','latex','FontSize',12);
        ylim([0,max(1,max(errorRe.spectral(i,:)))]); grid on; colororder(newcolors);
end
    leg = legend(txtleg2,'Orientation','Horizontal','FontSize',9);
    leg.Layout.Tile = 'north';
    sgtitle('Kappa Error in $\Omega(p)$ around '+txtbase,'Interpreter','latex','FontSize',16);

    %% Plot the percent error in the N^2 randomly drawn samples
fig5 = figure; tiledlayout(1,2);
nexttile
semilogy(real([errorRand.omega]),'*');
    title('$\Omega$ Relative Error','Interpreter','latex','FontSize',14);
    xlabel('Sample \#','Interpreter','latex','FontSize',14);
    ylabel('$\big|\frac{\Omega_{dielctric}-\Omega_{Xie}}{\Omega_{dielectric}}\big|$',...
        'Interpreter','latex','FontSize',14);
    xlim([0,N^2]); grid on
nexttile
semilogy(imag([errorRand.omega]),'*');
    title('$\gamma$ Relative Error','Interpreter','latex','FontSize',14);
    xlabel('Sample \#','Interpreter','latex','FontSize',14);
    ylabel('$\big|\frac{\gamma_{dielectric}-\gamma_{Xie}}{\gamma_{dielectric}}\big|$'...
        ,'Interpreter','latex','FontSize',14);
    xlim([0,N^2]); grid on

txt = ['Kappa ($\kappa$=',num2str(kappa),'), ',num2str(var*100),'\% variation on $(k,\sigma,\mu)$=(',...
    num2str(baseparams(1)),',',num2str(baseparams(2)),',',num2str(baseparams(3)),')'];
sgtitle(txt,'Interpreter','latex','FontSize',16)
    txt = ['Maximum Relative Error: ',num2str(errorMax,'%0.4e')];
    leg = legend('','Orientation','Horizontal','color','none'); legend boxoff;
    leg.Layout.Tile = 'north';
    annotation('textbox',[.25 .72 .6 .2], ...
        'String',txt,'EdgeColor','none','Interpreter','latex','FontSize',12)

%% Save the figures
% eps is a vector format file type that works well in LaTeX
% pdf is a vector format file type that is easy to share and view

savefig(fig1,  ['Figs\LevelCurvesKap' int2str(kappa) '_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.fig'])
hgexport(fig1, ['Figs\LevelCurvesKap' int2str(kappa) '_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
hgexport(fig1, ['Figs\LevelCurvesKap' int2str(kappa) '_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');
savefig(fig2,  ['Figs\LevelCurvesKap' int2str(kappa) '_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.fig'])
hgexport(fig2, ['Figs\LevelCurvesKap' int2str(kappa) '_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
hgexport(fig2, ['Figs\LevelCurvesKap' int2str(kappa) '_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');

savefig(fig3,  ['Figs\RelErrorKap' int2str(kappa) '_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.fig'])
hgexport(fig3, ['Figs\RelErrorKap' int2str(kappa) '_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
hgexport(fig3, ['Figs\RelErrorKap' int2str(kappa) '_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');
savefig(fig4,  ['Figs\RelErrorKap' int2str(kappa) '_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.fig'])
hgexport(fig4, ['Figs\RelErrorKap' int2str(kappa) '_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
hgexport(fig4, ['Figs\RelErrorKap' int2str(kappa) '_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');

savefig(fig5,  ['Figs\SampRandKap' int2str(kappa) '_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.fig'])
hgexport(fig5, ['Figs\SampRandKap' int2str(kappa) '_P' int2str(Nparams) '_N' int2str(N^2) '_var' num2str(var*100) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
hgexport(fig5, ['Figs\SampRandKap' int2str(kappa) '_P' int2str(Nparams) '_N' int2str(N^2) '_var' num2str(var*100) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');