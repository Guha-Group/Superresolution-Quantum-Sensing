% Determines the CFI of a marginalized BSPADE measurement where 
% marginalization occurs over a gaussian pointing error.

% sampling dimensionality for marginalization
edim = 101;

% rayleigh limit
sigma = 1;

% Total photon count
M = 1e4;

% separation range to consider
s = sigma.*exp(linspace(log(.01),log(1),101))';
s = permute(s,[3,2,1]);
ds = s(2)-s(1);

% container for CFI of the BSPADE measurement
%CFI_BS = zeros(numel(s)-1,M+1);
dm = 100;
M1_range = 1:dm:M;
M2_range = M-M1_range;
CFI_BS = zeros(numel(s),numel(M2_range));

%parpool(94)
%parfor n=1:numel(M1_range)
for n=1:numel(M1_range)
    M2 = M2_range(n);

    disp('Got to 1')
    
    % evaluate the distribution on the pointing error
    e_sig = sqrt(sigma^2/(M-M2));
    e = 5*e_sig*linspace(-1,1,edim)';
    e = permute(e,[2,1]);
    de = e(2)-e(1);
    p_e = normpdf(e,0,e_sig);


    % evaluate the conditional binomial distribution
    q0 = (0:M2)';
    Q = .5*(exp(-(e+s).^2/(2*sigma^2)) + exp(-(e-s).^2/(2*sigma^2)));
    disp('Got to 2')
    p_BS =  BSPADELikelihood(q0,M2,e,s,0,sigma,1);

    disp('Got to 3')
    % marginalize the binomial distribution with respect to pointing error
    P_BS = sum(p_BS.*p_e,2)*de;

    % compute the numerical derivative of the marginalized BS
    dP_BS = p_BS .* (q0-M2.*Q)./(1-Q).*(e.*tanh(e.*s/sigma^2)-s)/sigma^2;
    dP_BS = sum(dP_BS.*p_e,2)*de;
    disp('Got to 4')
    % compute the CFI (w analytic gradient)
    CFI_BS(:,n) = squeeze(sum( dP_BS.^2 ./ (P_BS+1e-20),1));

    % display count
    fprintf('%d : %d',n,M);

end

save('Calibration_OptimalSwitching_FisherInformation.mat');
%,'CFI_BS','M','M2_range','sigma','s','dm','edim')

%{
CFI_BS_opt = max(CFI_BS,[],2);


% CFI of direct imaging
x = 5*sigma*linspace(-1,1,1001);
dx = x(2)-x(1);
p0 = @(y) normpdf(x-y,0,sigma);
dp0 = @(y) -(x-y).*normpdf(x-y,0,sigma)/sigma^2;
p_DD = (p0(-s) + p0(s))/2;
dp_DD = (dp0(-s)-dp0(s))/2;
CFI_DD = sum(dp_DD.^2./(p_DD+1e-20))*dx;

s = squeeze(s);
CFI_DD = squeeze(CFI_DD);


% total CFI
M1 = M-M2_range;
CFI_tot = M1.*CFI_DD + CFI_BS;
CRB = CFI_tot.^(-1);


% optimal photon alocation between DD and BSPADE during Calibration Stage
[CRB_opt,M2_opt] = min(CRB,[],2);
M1_opt = M-M2_opt;

figure
t=tiledlayout(1,2);

nexttile(1)
plot(s/sigma,M1_opt/M,'k','LineWidth',2)
set(gca,'xscale','log')
ylabel('Optimal DI Allocation $M_1/M$','interpreter','latex')
xlabel('$s/\sigma$','interpreter','latex')
axis square
box on

nexttile(2)
plot(s/sigma,CRB_opt/sigma^2,'k','LineWidth',2);
set(gca,'xscale','log')
ylabel('CRB','interpreter','latex')
xlabel('$s/\sigma$','interpreter','latex')
title(t,{'Optimal Switching',sprintf('$M=%d$',M)},'interpreter','latex')
axis square
box on

% plot the CFI for both measurements
figure
t=tiledlayout(1,3);

nexttile(1)
hold on
plot(s/sigma,CFI_DD,'--k','LineWidth',2)
plot(s/sigma,CFI_BS_opt,'k','LineWidth',2)
hold off
xlabel('$s/\sigma$','interpreter','latex')
ylabel('Direct Imaging CFI $\mathcal{I}_{DD}^{(s)}$','interpreter','latex')
set(gca,'xscale','log')
axis square
box on


nexttile(2)
colormap('winter')
[S,M1_ratio] = meshgrid(s,flip(M1/M));
surf(S,M1_ratio,CFI_BS)
shading interp
set(gca,'xscale','log')
xlabel('$s/\sigma$','interpreter','latex')
ylabel('Splitting Ratio $M_1/M$','interpreter','latex')
zlabel('CFI BSPADE $\mathcal{I}_{BS}^{(s)}$','interpreter','latex')
axis square

nexttile(3)
contour(S,M1_ratio,CFI_BS,15)
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
%}

