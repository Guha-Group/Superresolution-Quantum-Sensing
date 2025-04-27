function fig = PlotReceiver(outputdata,receiver)
    % A function for plotting receiver metadata to interpret performance.
    % The receiver list could be extended if desired 

    validReceivers = {'DirectImaging','StaticBSPADE','AdaptiveBSPADE','ThresholdBSPADE'};
    validateReceiver = @(x) any(validatestring(x,validReceivers));

    % default inputs
    default_sigma = 1; % sets the resolution scale

    % input parser
    p = inputParser;
    addRequired(p,'outputdata')
    addRequired(p,'receiver',validateReceiver)
    parse(p,outputdata,receiver)
    
    switch receiver
        case 'DirectImaging'
            fig = PlotDirectImaging(outputdata);
        case 'StaticBSPADE'
            fig = PlotStaticBSPADE(outputdata);
        case 'AdaptiveBSPADE'
            fig = PlotAdaptiveBSPADE(outputdata);
        case 'ThresholdBSPADE'
            fig = PlotThresholdBSPADE(outputdata);
    end
end

%% PLOTTING FUNCTIONS
function fig = PlotDirectImaging(outputdata)
    sigma = 1;
    % extract structure params
    s = outputdata.s;
    k = outputdata.k;
    xsk_0 = outputdata.xsk_0;
    photons = outputdata.photons;
    xsk_est = outputdata.xsk_est;
    P_s = outputdata.P_s;
    P_k = outputdata.P_k;

    % extract particular values
    s0 = xsk_0(2);      k0 = xsk_0(3);
    s_est = xsk_est(2); k_est = xsk_est(3);
    M = photons(1);
    N = photons(3);

    % some plotting colors
    c1 = [0 0.4470 0.7410];
    c2 = [.1 0.7410 0.4470];
    c3 = [0.8410 0.2 .1];

    % some objects for shadding below functions (just visual snazziness) 
    ss = [squeeze(s); flip(squeeze(s))];
    kk = [squeeze(k); flip(squeeze(k))];
    f2 = [zeros(numel(s),1); flip(squeeze(P_s))];
    f3 = [zeros(numel(k),1); flip(squeeze(P_k))];
    
    % plot
    fig = figure;
    t=tiledlayout(1,2,"TileSpacing","compact","Padding","compact");
    
    % plot the distribution on the separation estimate for direct imaging
    nexttile(1)
    hold on
    plot(squeeze(s)/sigma,squeeze(P_s),'Color',c1,'LineWidth',1.5)
    xline(s_est,'k','LineWidth',1.5)
    xline(s0,'--k','LineWidth',1.5)
    fill(ss,f2,c1,'FaceAlpha',.1,'EdgeColor',c1)
    hold off
    axis square
    box on
    xlabel('$s/\sigma$','interpreter','latex')
    ylabel({'Separation Probability Density','$p(s)$'},'interpreter','latex')
    legend({'Prior $\tilde{p}(s|\mathbf{x})$ (Direct Imaging)',...
        ['MMSE Estimate $\check{s}/\sigma=',sprintf('%.3f$',s_est/sigma)],...
        ['Ground Truth $s/\sigma=',sprintf('%.3f$',s0)]},'interpreter','latex')
    title({'\textbf{CALIBRATION STAGE}',...
        'Separation Estimation',...
        'Photon Allocation:',sprintf('$M=%d$',M)},'interpreter','latex')
    
    % plot the distribution on the brightness estimate for direct imaging
    nexttile(2)
    hold on
    plot(squeeze(k),squeeze(P_k),'Color',c3,'LineWidth',1.5)
    xline(k_est,'k','LineWidth',1.5)
    xline(k0,'--k','LineWidth',1.5)
    fill(kk, f3,'r','FaceAlpha',.1,'EdgeColor',c3)
    hold off
    axis square
    box on
    xlabel('$\kappa$','interpreter','latex')
    ylabel('Brightness Probability Density $p(k)$','interpreter','latex')
    legend({'Marginal Posterior $p(\kappa|\mathbf{x}'',q,\mathbf{x})$',...
        ['MMSE Estimate $\check{\kappa}=',sprintf('%.2f$',k_est)],...
        ['Ground Truth $\kappa=',sprintf('%.2f$',k0)]},'interpreter','latex')
    title({'\textbf{SENSING STAGE}',...
           'Brightness Estimation',...
           'Photon Allocation:',sprintf('$N=%d$',N)},'interpreter','latex')
    title(t,'DirectImaging','interpreter','latex')
end

function fig = PlotStaticBSPADE(outputdata)
    sigma = 1;
    % extract structure params
    s = outputdata.s;
    k = outputdata.k;
    xsk_0 = outputdata.xsk_0;
    photons = outputdata.photons;
    xsk_est = outputdata.xsk_est;
    p_s = outputdata.p_s;
    P_s = outputdata.P_s;
    P_k = outputdata.P_k;

    % extract particular values
    s0 = xsk_0(2);      k0 = xsk_0(3);
    s_est = xsk_est(2); k_est = xsk_est(3);
    M1 = photons(1);
    M2 = photons(2);
    N = photons(3);

    % some plotting colors
    c1 = [0 0.4470 0.7410];
    c2 = [.1 0.7410 0.4470];
    c3 = [0.8410 0.2 .1];

    % some objects for shadding below functions (just visual snazziness) 
    ss = [squeeze(s); flip(squeeze(s))];
    kk = [squeeze(k); flip(squeeze(k))];
    f1 = [zeros(numel(s),1); flip(squeeze(p_s))];
    f2 = [zeros(numel(s),1); flip(squeeze(P_s))];
    f3 = [zeros(numel(k),1); flip(squeeze(P_k))];
    
    % plot
    fig = figure;
    t = tiledlayout(1,2,"TileSpacing","compact","Padding","compact");
    
    % plot posterior update before and after BSPADE measurement
    nexttile(1)
    hold on
    plot(squeeze(s)/sigma,squeeze(p_s),'Color',c1,'LineWidth',1.5)
    plot(squeeze(s)/sigma,squeeze(P_s),'Color',c2,'LineWidth',1.5)
    xline(s_est,'k','LineWidth',1.5)
    xline(s0,'--k','LineWidth',1.5)
    fill(ss,f1,c1,'FaceAlpha',.1,'EdgeColor',c1)
    fill(ss,f2,c2,'FaceAlpha',.1,'EdgeColor',c2)
    hold off
    axis square
    box on
    xlabel('$s/\sigma$','interpreter','latex')
    ylabel({'Separation Probability Density','$p(s)$'},'interpreter','latex')
    legend({'Prior $\tilde{p}(s|\mathbf{x})$ (Direct Imaging)',...
        'Marginal Posterior $p(s|q,\mathbf{x})$ (BSPADE)',...
        ['MMSE Estimate $\check{s}/\sigma=',sprintf('%.3f$',s_est/sigma)],...
        ['Ground Truth $s/\sigma=',sprintf('%.3f$',s0)]},'interpreter','latex')
    title({'\textbf{CALIBRATION STAGE}',...
        'Separation Estimation',...
        'Photon Allocation:',sprintf('$[M_{1},M_{2}]=[%d,%d]$',M1,M2)},'interpreter','latex')
    
    
    % plot the posterior on the brightness parameter
    nexttile(2)
    hold on
    plot(squeeze(k),squeeze(P_k),'Color',c3,'LineWidth',1.5)
    xline(k_est,'k','LineWidth',1.5)
    xline(k0,'--k','LineWidth',1.5)
    fill(kk, f3,'r','FaceAlpha',.1,'EdgeColor',c3)
    hold off
    axis square
    box on
    xlabel('$\kappa$','interpreter','latex')
    ylabel('Brightness Probability Density $p(k)$','interpreter','latex')
    legend({'Marginal Posterior $p(\kappa|\mathbf{x}'',q,\mathbf{x})$',...
        ['MMSE Estimate $\check{\kappa}=',sprintf('%.2f$',k_est)],...
        ['Ground Truth $\kappa=',sprintf('%.2f$',k0)]},'interpreter','latex')
    title({'\textbf{SENSING STAGE}',...
           'Brightness Estimation',...
           'Photon Allocation:',sprintf('$N=%d$',N)},'interpreter','latex')
    
end
function fig = PlotAdaptiveBSPADE(outputdata)
    sigma = 1;
    % extract structure params
    s = outputdata.s;
    k = outputdata.k;
    xsk_0 = outputdata.xsk_0;
    photons = outputdata.photons;
    xsk_est = outputdata.xsk_est;
    p_s = outputdata.p_s;
    P_s = outputdata.P_s;
    P_k = outputdata.P_k;
    dm = outputdata.dm;
    iter = outputdata.iter;
    expected_var_s = outputdata.expected_var_s;

    % extract particular values
    s0 = xsk_0(2);      k0 = xsk_0(3);
    s_est = xsk_est(2); k_est = xsk_est(3);
    M1 = photons(1);
    M2 = photons(2);
    N = photons(3);


    % some plotting colors
    c1 = [0 0.4470 0.7410];
    c2 = [.1 0.7410 0.4470];
    c3 = [0.8410 0.2 .1];

    % some objects for shadding below functions (just visual snazziness) 
    ss = [squeeze(s); flip(squeeze(s))];
    kk = [squeeze(k); flip(squeeze(k))];
    f1 = [zeros(numel(s),1); flip(squeeze(p_s))];
    f2 = [zeros(numel(s),1); flip(squeeze(P_s))];
    f3 = [zeros(numel(k),1); flip(squeeze(P_k))];
    
    % plot
    fig = figure;
    t = tiledlayout(1,3,"TileSpacing","compact","Padding","compact");
    
    % plot posterior update before and after BSPADE measurement
    nexttile(1)
    hold on
    plot(squeeze(s)/sigma,squeeze(p_s),'Color',c1,'LineWidth',1.5)
    plot(squeeze(s)/sigma,squeeze(P_s),'Color',c2,'LineWidth',1.5)
    xline(s_est,'k','LineWidth',1.5)
    xline(s0,'--k','LineWidth',1.5)
    fill(ss,f1,c1,'FaceAlpha',.1,'EdgeColor',c1)
    fill(ss,f2,c2,'FaceAlpha',.1,'EdgeColor',c2)
    hold off
    axis square
    box on
    xlabel('$s/\sigma$','interpreter','latex')
    ylabel({'Separation Probability Density','$p(s)$'},'interpreter','latex')
    legend({'Prior $\tilde{p}(s|\mathbf{x})$ (Direct Imaging)',...
        'Marginal Posterior $p(s|q,\mathbf{x})$ (BSPADE)',...
        ['MMSE Estimate $\check{s}/\sigma=',sprintf('%.3f$',s_est/sigma)],...
        ['Ground Truth $s/\sigma=',sprintf('%.3f$',s0)]},'interpreter','latex')
    title({'\textbf{CALIBRATION STAGE}',...
        'Separation Estimation',...
        'Photon Allocation:',sprintf('$[M_{1},M_{2}]=[%d,%d]$',M1,M2)},'interpreter','latex')
    
    
    % plot the posterior on the brightness parameter
    nexttile(2)
    hold on
    plot(squeeze(k),squeeze(P_k),'Color',c3,'LineWidth',1.5)
    xline(k_est,'k','LineWidth',1.5)
    xline(k0,'--k','LineWidth',1.5)
    fill(kk, f3,'r','FaceAlpha',.1,'EdgeColor',c3)
    hold off
    axis square
    box on
    xlabel('$\kappa$','interpreter','latex')
    ylabel('Brightness Probability Density $p(k)$','interpreter','latex')
    legend({'Marginal Posterior $p(\kappa|\mathbf{x}'',q,\mathbf{x})$',...
        ['MMSE Estimate $\check{\kappa}=',sprintf('%.2f$',k_est)],...
        ['Ground Truth $\kappa=',sprintf('%.2f$',k0)]},'interpreter','latex')
    title({'\textbf{SENSING STAGE}',...
           'Brightness Estimation',...
           'Photon Allocation:',sprintf('$N=%d$',N)},'interpreter','latex')

    
    % plot switching fraction instance
    nexttile(3)
    hold on
    scatter(dm:dm:iter*dm,expected_var_s(1:iter),'k','filled')
    plot(dm:dm:iter*dm,expected_var_s(1:iter),'k','LineWidth',1.5)
    hold off
    axis square
    box on
    ylabel('Expected Separation Variance','interpreter','latex')
    xlabel('Number of Photons $M_1$','interpreter','latex');
    xlim([dm,M])
    set(gca,'xscale','log')
    title({'\textbf{ADAPTIVE SWITCHING}','Trace of Dynamic Cost Function','',''},'interpreter','latex')
end

function fig = PlotThresholdBSPADE(metadata)
    fig = PlotStaticBSPADE(metadata);
end
