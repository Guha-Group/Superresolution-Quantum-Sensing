%% FIGURES
figure
if adaptive
    tiledlayout(1,3,"TileSpacing","compact","Padding","compact")
else
    tiledlayout(1,2,"TileSpacing","compact","Padding","compact")
end

ss = [squeeze(s);fliplr(squeeze(s))];
kk = [squeeze(k); flip(squeeze(k))];
f1 = [zeros(numel(s),1);squeeze(p_s)];
f2 = [zeros(numel(s),1);squeeze(P_s)];tiledlayout(1,3,"TileSpacing","compact","Padding","compact")
f3 = [zeros(numel(k),1);flip(squeeze(P_k))];

% plot posterior update before and after BSPADE measurement
nexttile(1)
hold on
plot(squeeze(s)/sigma,squeeze(p_s),'Color',[0 0.4470 0.7410],'LineWidth',1.5)
plot(squeeze(s)/sigma,squeeze(P_s),'Color',[.1 0.7410 0.4470 ],'LineWidth',1.5)
xline(s_mmse,'k','LineWidth',1.5)
xline(s0,'--k','LineWidth',1.5)
fill(ss,f1,[0 0.4470 0.7410],'FaceAlpha',.1,'EdgeColor',[0 0.4470 0.7410])
fill(ss,f2,[.1 0.7410 0.4470 ],'FaceAlpha',.1,'EdgeColor',[.1 0.7410 0.4470 ])
hold off
axis square
box on
xlabel('$s/\sigma$','interpreter','latex')
ylabel({'Separation Probability Density','$p(s)$'},'interpreter','latex')
legend({'Prior $\tilde{p}(s|\mathbf{x})$ (Direct Imaging)',...
    'Marginal Posterior $p(s|q,\mathbf{x})$ (BSPADE)',...
    ['MMSE Estimate $\check{s}/\sigma=',sprintf('%.3f$',s_mmse/sigma)],...
    ['Ground Truth $s/\sigma=',sprintf('%.3f$',s0)]},'interpreter','latex')
title({'\textbf{CALIBRATION STAGE}',...
    'Separation Estimation',...
    'Photon Allocation:',sprintf('$[M_{1},M_{2}]=[%d,%d]$',M1,M2)},'interpreter','latex')


% plot the posterior on the brightness parameter
nexttile(2)
hold on
plot(squeeze(k),squeeze(P_k),'Color',[0.8410 0.2 .1],'LineWidth',1.5)
xline(k_mmse,'k','LineWidth',1.5)
xline(k0,'--k','LineWidth',1.5)
fill(kk, f3,'r','FaceAlpha',.1,'EdgeColor',[0.8410 0.2 .1])
hold off
axis square
box on
xlabel('$\kappa$','interpreter','latex')
ylabel('Brightness Probability Density $p(k)$','interpreter','latex')
legend({'Marginal Posterior $p(\kappa|\mathbf{x}'',q,\mathbf{x})$',...
    ['MMSE Estimate $\check{\kappa}=',sprintf('%.2f$',k_mmse)],...
    ['Ground Truth $\kappa=',sprintf('%.2f$',k0)]},'interpreter','latex')
title({'\textbf{SENSING STAGE}',...
       'Brightness Estimation',...
       'Photon Allocation:',sprintf('$N=%d$',N)},'interpreter','latex')


% plot switching fraction instance
if adaptive
    nexttile(3)
    hold on
    scatter(dm:dm:iter*dm,expected_var_s,'k','filled')
    plot(dm:dm:iter*dm,expected_var_s,'k','LineWidth',1.5)
    hold off
    axis square
    box on
    ylabel('Expected Separation Variance','interpreter','latex')
    xlabel('Number of Photons $M_1$','interpreter','latex');
    xlim([dm,M])
    set(gca,'xscale','log')
    title({'\textbf{ADAPTIVE SWITCHING}','Trace of Dynamic Cost Function','',''},'interpreter','latex')
end
