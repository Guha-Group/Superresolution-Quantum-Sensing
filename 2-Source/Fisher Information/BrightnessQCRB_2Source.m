%% QFIM for x0,s,kappa state
sigma = 1;
s = 10.^(linspace(-3,1,1000))/sigma;
phi = exp(-1/2*(s/sigma).^2);
g = (s.*phi/sigma).^2;
kappa = (0:.1:.4)';
tau= 4*kappa.^2;
coeffvct = [tau.^2, zeros(size(tau)),-2*tau,-2*(1-tau),ones(size(tau))];

N = 1;
QCRB_kappa = (1/4)*(1-4*kappa.^2).*(g-1)...
    ./(g+phi.^2 -1 )/N;

plot(s/sigma,QCRB_kappa,'LineWidth',1.5);
set(gca,'yscale','log')
set(gca,'xscale','log')

xlabel('Source Separation $s/\sigma$','interpreter','latex')
ylabel('Brightness QCRB$(\check{\kappa})$','interpreter','latex')
legend_names = arrayfun(@(j) sprintf('%.2f',kappa(j)),1:numel(kappa),'UniformOutput',false);
leg = legend(legend_names);


%% QFIM for cumulative Estimation and Sensing Stages

sigma = 1;
s = 10.^(linspace(-3,1,1000))/sigma;
phi = exp(-1/2*(s/sigma).^2);
gamma = (s.*phi/sigma).^2;
kappa = (0:.1:.4)';

alpha=zeros(size(kappa));
for k = 1:numel(kappa)
    rts = roots(coeffvct(k,:));
    alpha(k) = rts(0<rts & rts<1);
end
tau = 4*kappa.^2;

numerator = (1-tau).*((1-gamma)-tau.*alpha.*(alpha-gamma));
denominator = 4.*alpha.*((1-tau.*alpha.^2).*(1-phi.^2)-gamma.*(1-tau.*alpha-(1-alpha).*phi.^2));
QCRB_kappa = numerator./denominator;
%QCRB_kappa_subrl_approx = (sigma./s).^2.*((1-tau.*alpha^2).*(1-tau))./(4*alpha.*((1-tau.*alpha.^2) - alpha.*(1-tau)));
QCRB_kappa_subrl_approx = (sigma./s).^2.*(1./(4*alpha.^2)).*((1./(alpha.*(1-tau)) - 1./(1-tau.*alpha.^2)).^(-1));

figure
hold on
plot(s/sigma,QCRB_kappa,'LineWidth',1.5);
%plot(s/sigma,QCRB_kappa_subrl_approx,'--','LineWidth',1.5);
hold off
set(gca,'yscale','log')
set(gca,'xscale','log')
legend
legend_names = arrayfun(@(j) sprintf('%.2f',kappa(j)),1:numel(kappa),'UniformOutput',false);
leg = legend(legend_names);
title({'Quantum Cramer-Rao Bound for Brightness Estimation','@ Optimal Splitting $\beta^*$'},'interpreter','latex')
title(leg,'$\kappa$','interpreter','latex')
xlabel('Source Separation $s/\sigma$','interpreter','latex')
ylabel({'$(N+M) \cdot QCRB(\check{\kappa}) $'},'interpreter','latex')



% plot the optimal allocation ratio as a function of separation
kappa = linspace(0+1e-5,.5-1e-5,10000)';
tau= 4*kappa.^2;
coeffvct = [tau.^2, zeros(size(tau)),-2*tau,-2*(1-tau),ones(size(tau))];
alpha_opt=zeros(size(kappa));
for k = 1:numel(kappa)
    rts = roots(coeffvct(k,:));
    alpha_opt(k) = rts(0<rts & rts<1);
end

figure
plot(kappa,alpha_opt,'k','LineWidth',1.5)
axis square
xlim([0,0.5]);
ylim([0,1])
xlabel('$|\kappa|$','interpreter','latex')
ylabel('$\beta^{*}$','interpreter','latex')
title('Optimal Calibration Stage Allocation Ratio','interpreter','latex')
