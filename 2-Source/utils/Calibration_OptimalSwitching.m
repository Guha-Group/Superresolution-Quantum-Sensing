% Determines the CFI of a marginalized BSPADE measurement where 
% marginalization occurs over a gaussian pointing error.
%{
% sampling dimensionality for marginalization
edim = 101;

% rayleigh limit
sigma = 1;

% Total photon count
M = 1e3;

% separation range to consider
s = sigma.*exp(linspace(log(.01),log(10),101))';
s = permute(s,[3,2,1]);
ds = s(2)-s(1);

% container for CFI of the BSPADE measurement
%CFI_BS = zeros(numel(s)-1,M+1);
dm = 1;
M1_range = 1:dm:M;
M2_range = M-M1_range;
CFI_BS = zeros(numel(s),numel(M2_range));

%parpool(94)
%parfor n=1:numel(M1_range)
for n=1:numel(M1_range)
    M1 = M1_range(n);
    M2 = M2_range(n);

    % evaluate the distribution on the pointing error
    e_sig = sqrt(sigma^2/M1);
    e = 5*e_sig*linspace(-1,1,edim)';
    e = permute(e,[2,1]);
    de = e(2)-e(1);
    p_e = normpdf(e,0,e_sig);

    % evaluate the conditional binomial distribution
    q0 = (0:dm:M2)';
    Q = exp(-(e.^2+s.^2)./(4*sigma^2)).*cosh(e.*s./(4*sigma^2));
    p_BS =  BSPADELikelihood(q0,M2,e,s,0,sigma,1);

    % marginalize the binomial distribution with respect to pointing error
    P_BS = sum(p_BS.*p_e,2)*de;

    % compute the derivative of the marginalized BS
    dP_BS = p_BS .* (q0-M2.*Q)./(1-Q).*(e.*tanh(e.*s/(4*sigma^2))-2*s)/(4*sigma^2);
    dP_BS = sum(dP_BS.*p_e,2)*de;

    % compute the CFI (w analytic gradient)
    CFI_BS(:,n) = squeeze(sum( dP_BS.^2 ./ (P_BS+1e-20),1))*dm;

    % display count
    fprintf('%d : %d\n',n,numel(M1_range));

end

save('Calibration_OptimalSwitching_FisherInformation.mat');
%,'CFI_BS','M','M2_range','sigma','s','dm','edim')
%}
load('Calibration_OptimalSwitching_FisherInformation.mat')

% CFI of direct imaging
x = 25*sigma*linspace(-1,1,1001);
dx = x(2)-x(1);
p0 = @(y) normpdf(x-y,0,sigma);
dp0 = @(y) -(x-y).*normpdf(x-y,0,sigma)/sigma^2;
p_DD = (p0(-s) + p0(s))/2;
dp_DD = (dp0(-s)-dp0(s))/2;
CFI_DD_per_photon = sum(dp_DD.^2./(p_DD+1e-20))*dx;

% return dimensions of objects
s = squeeze(s);
CFI_DD_per_photon = squeeze(CFI_DD_per_photon);

% total CFI
CFI_DD = M1_range.*CFI_DD_per_photon;
CFI_tot = CFI_DD + CFI_BS;
CRB = 1./CFI_tot;

% optimal photon alocation between DD and BSPADE during Calibration Stage
[CFI_opt,opt_id] = max(CFI_tot,[],2);
M1_opt = M1_range(opt_id)';
CRB_opt = 1./CFI_opt;

% get optimal CFI contributions from direct imaging and BSPADE
CFI_BS_opt = CFI_BS(sub2ind(size(CFI_BS),(1:numel(s))',opt_id));
CFI_DD_opt = M1_opt.*CFI_DD_per_photon;

%% FIGURES
figure
t=tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

nexttile(1)
plot(s/sigma,M1_opt/M,'k','LineWidth',2)
set(gca,'xscale','log')
ylabel('Optimal DI Allocation $M_1/M$','interpreter','latex')
xlabel('$s/\sigma$','interpreter','latex')
axis square
box on

nexttile(2)
hold on
plot(s/sigma,ones(size(s))/M,'--k','LineWidth',1.5)
plot(s/sigma,CRB_opt/sigma^2,'k','LineWidth',1.5);
plot(s/sigma,1./(M*CFI_DD_per_photon),'--r','LineWidth',1.5)
hold off
set(gca,'xscale','log')
set(gca,'yscale','log')
ylabel('CRB','interpreter','latex')
xlabel('$s/\sigma$','interpreter','latex')
legend({'2-Stage CRB','QCRB','Direct Imaging'},'interpreter','latex')
title(t,{'Optimal Switching',sprintf('$M=%d$',M)},'interpreter','latex')
axis square
box on

nexttile(3)
hold on
plot(s/sigma,M*ones(size(s)),'--k','LineWidth',1.5)
plot(s/sigma,CFI_opt,'k','LineWidth',1.5)
plot(s/sigma,CFI_BS_opt,'--b','LineWidth',1.5)
plot(s/sigma,CFI_DD_opt,'--r','LineWidth',1.5)
hold off
xlabel('$s/\sigma$','interpreter','latex')
ylabel('CFI $\mathcal{I}^{(s)}$','interpreter','latex')
set(gca,'xscale','log')
set(gca,'yscale','log')
legend({'QFI','2-Stage Optimal','BSPADE Contribution','Direct Imaging Contribution'},'interpreter','latex')
title('Fisher Information','interpreter','latex')
axis square
box on

% plot the change in CFI or CRB as a function of separation and M1
% allocation for each component
figure
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')

% plot CFI chart of Direct Imaging
nexttile(1)
colormap(turbo)
imagesc([min(s),max(s)],[min(M1_range),max(M1_range)],CFI_DD');
axis square
set(gca,'ydir','normal')
xlabel('$s/\sigma$','interpreter','latex')
ylabel('$M_1$','interpreter','latex')
cbar=colorbar;
ylabel(cbar,'$M_1 CFI^{DD}(s)$','interpreter','latex')
clim([min(CFI_tot(:)),max(CFI_tot(:))])
set(gca,'ColorScale','log')
title('Direct Imaging')


% plot CFI chart of BSPADE
nexttile(2)
imagesc([min(s),max(s)],[min(M1_range),max(M1_range)],CFI_BS');
axis square
set(gca,'ydir','normal')
xlabel('$s/\sigma$','interpreter','latex')
ylabel('$M_1$','interpreter','latex')
cbar=colorbar;
ylabel(cbar,'$CFI^{2-Stage}_{s}(s,M_1)$','interpreter','latex')
clim([min(CFI_tot(:)),max(CFI_tot(:))])
set(gca,'ColorScale','log')
title('BSPADE')

% plot cumulative CFI chart
nexttile(3)
imagesc([min(s),max(s)],[min(M1_range),max(M1_range)],CFI_tot');
axis square
set(gca,'ydir','normal')
xlabel('$s/\sigma$','interpreter','latex')
ylabel('$M_1$','interpreter','latex')
cbar=colorbar;
ylabel(cbar,'$CFI^{tot}_{s}(s,M_1) = M_1 CFI^{DD}(s) + CFI^{2-Stage}(s,M_1)$','interpreter','latex')
clim([min(CFI_tot(:)),max(CFI_tot(:))])
set(gca,'ColorScale','log')
title('Cumulative')


% plot cumulative CRB chart
nexttile(4)
imagesc([min(s),max(s)],[min(M1_range),max(M1_range)],CRB');
axis square
set(gca,'ydir','normal')
xlabel('$s/\sigma$','interpreter','latex')
ylabel('$M_1$','interpreter','latex')
cbar=colorbar;
ylabel(cbar,'$1/CFI^{tot}_{s}(s,M_1)$','interpreter','latex')
set(gca,'ColorScale','log')
title('CRB')


%{
% plot the CFI for both measurements
figure
t=tiledlayout(1,3);

nexttile(2)
colormap('winter')
[S,M1_ratio] = meshgrid(s,M1_range/M);
surf(S,M1_ratio,CFI_BS')
shading interp
set(gca,'xscale','log')
xlabel('$s/\sigma$','interpreter','latex')
ylabel('Splitting Ratio $M_1/M$','interpreter','latex')
zlabel('CFI BSPADE $\mathcal{I}_{BS}^{(s)}$','interpreter','latex')
axis square

nexttile(3)
contour(S,M1_ratio,CFI_BS',15)
axis square
set(gca,'xscale','log')



% Make a contour plot of the log CFI
figure
hold on
imagesc([min(s),max(s)],[min(M1_range),max(M1_range)],log(CFI_tot'),'AlphaData',.25)
contour(s,M1_range,log(CFI_tot'),30)
hold off
set(gca,'ydir','normal')
axis square
xlabel('$s/\sigma$','Interpreter','latex')
ylabel('$M_1$','Interpreter','latex')


% Test BSPADE CFIM with perfect pointing
q0 = (0:M)';
s = sigma.*exp(linspace(log(.01),log(10),101))';
s = permute(s,[3,2,1]);
p_BS = BSPADELikelihood(q0,M,0,s,0,sigma,0);
Q = exp(-(s/2*sigma).^2);
% promote array dimensions
%Q0 = repmat(q0,size(Q)); P0 = repmat(Q,size(q0));
%p_BS = binopdf(Q0,M,P0);
CFI_BS_ref = sum(p_BS.*((q0-M.*Q)./(1-Q).*(-s/2/sigma^2)).^2,1);
CFI_BS_analytic = M*(s/2/sigma^2).^2.*(Q)./(1-Q);
figure
hold on
plot(squeeze(s),squeeze(CFI_BS_ref),'r')
plot(squeeze(s),M*ones(1,numel(s)),'k')
plot(squeeze(s),squeeze(CFI_BS_analytic),'--b')
set(gca,'xscale','log')
set(gca,'yscale','log')
hold off
%}