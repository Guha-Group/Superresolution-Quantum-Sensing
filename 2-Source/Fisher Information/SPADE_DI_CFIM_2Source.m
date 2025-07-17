%% A script for computing the Quantum Fisher information matrix for the
% two source problem as well as the classical Fisher information for 
% direct imaging and HG-SPADE.
% 
% Author(s): Nico Deshler

% imaging system
sigma = 1;      % rayleigh limit for gaussian PSF

% parameters (as column vectors)
x0 = sigma*(0:0.02:0.2)';            % geometric midpoint   
s = sigma*10.^linspace(-3,0,100)';   % half-separation
kappa = (0:0.05:.45)';               % brightness bias

% assign array dimension to each parameter
x0 = permute(x0,[2,1]);             % dim 2
s = permute(s,[3,2,1]);             % dim 3
kappa = permute(kappa,[4,3,2,1]);   % dim 4

% transformed parameters
b1 = 1/2 - kappa;
b2 = 1/2 + kappa;
x1 = x0-s; 
x2 = x0+s;

% list of figures to plot
plot_list = [1,1,0,0,1,0,0,0,0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% QUANTUM BOUND %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the QFIM
Q = zeros(3,3,1,numel(s),numel(kappa));
phi = exp(-(s/sigma).^2/2);
Q(1,1,:,:,:) = 1 - (phi.*s/sigma).^2.*(1-4*kappa.^2);
Q(1,2,:,:,:) = 2.*kappa.*ones(size(s));
Q(1,3,:,:,:) = 4*s.*phi.^2.*ones(size(kappa));
Q(2,2,:,:,:) = 1;
Q(3,3,:,:,:) = 4*sigma^2.*(1-phi.^2)./(1-4.*kappa.^2);
Q = Q + pagetranspose(Q.*(1-eye(3)));
Q(1:2,1:2,:,:,:) = Q(1:2,1:2,:,:,:)/sigma;

% calculate quantum cramer-rao matrix
QCRB = pageinv(Q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% DIRECT IMAGING %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the direct imaging CFIM
% imaging detector plane
x = sigma*linspace(-20,20,10001)';
dx = x(2)-x(1);

% direct imaging pdf
p0 =  @(u) normpdf(u,0,sigma);
dp0 = @(u) -(u/sigma^2).*normpdf(u,0,sigma);
px = b1.*p0(x-x1)+ b2.*p0(x-x2);

% partial derivatives with respect to each parameter 
dpx_1 = -(b1.*dp0(x-x1) + b2.*dp0(x-x2));                       % dp_dx0
dpx_2 = b1.*dp0(x-x1) - b2.*dp0(x-x2);                          % dp_ds
dpx_3 = (-p0(x-x1) + p0(x-x2)).*ones(size(kappa));              % dp_dk

% collect all partial derivatives in an array
Dpx = cat(5,dpx_1,dpx_2,dpx_3);

% compute the entries of the CFIM for direct imaging
I_DD = zeros([3,3,numel(x0),numel(s),numel(kappa)]); % psf is shift invariant so IDD won't ever depend on x0
for i = 1:3
    for j = i:3
        I_DD(i,j,:,:,:) = sum(Dpx(:,:,:,:,i).*Dpx(:,:,:,:,j)./px,1)*dx;
    end
end
I_DD = I_DD + pagetranspose(I_DD.*(1-eye(3)));
I_DD(1:2,1:2,:,:,:) = I_DD(1:2,1:2,:,:,:)*sigma^2; % normalize physical dimensions

% compute direct imaging cramer-rao matrix
CRB_DD = pageinv(I_DD);

% MATHEMATICA OUTPUT
% compute the CFIM in the sub-diffraction taylor approximation
I_DD_approx = zeros(3,3,1,numel(s),numel(kappa));
I_DD_approx(1,1,1,:,:) = (1-(s/sigma).^2 .*(1-4*kappa.^2) +((s/sigma).^2 .*(1-4*kappa.^2)).^2 -(s/sigma).^6 .*(1-4*kappa.^2).^2.*(1-12*kappa.^2));
I_DD_approx(1,2,1,:,:) = 2*kappa.*(1-(s/sigma).^2 .*(1-4*kappa.^2));
I_DD_approx(1,3,1,:,:) = 2*s.*(1-(s/sigma).^2 .*(1-4*kappa.^2));
I_DD_approx(2,2,1,:,:) = (4*kappa.^2 + 2*(s/sigma).^2.*(1-10*kappa.^2+24.*kappa.^4));
I_DD_approx(2,3,1,:,:) = 4*s.*kappa.*(1-3*(s/sigma).^2.*(1-4*kappa.^2));
I_DD_approx(3,3,1,:,:) = 4*s.^2.*ones(size(kappa));
I_DD_approx = I_DD_approx/(sigma^2);
I_DD_approx = repmat(I_DD_approx,[1,1,numel(x0),1,1]);
I_DD_approx = I_DD_approx + pagetranspose(I_DD_approx.*(1-eye(3)));
I_DD_approx(1:2,1:2,:,:,:) = I_DD_approx(1:2,1:2,:,:,:)*sigma^2; % normalize physical dimensions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% HG-SPADE %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% total mode count
q = (0:1e3)';

% transformed parameters
b1 = 1/2 - kappa;
b2 = 1/2 + kappa;
x1 = x0-s; 
x2 = x0+s;

% HG-SPADE PMF
pq0 = @(u) exp(-(u/2/sigma).^2) .* (u/2/sigma).^(2*q)./ factorial(q); 
dpq0 = @(u) (q.*(2*sigma./u).^2 - 1).*(u/2/sigma^2).*pq0(u);
pq = b1.*pq0(x1) + b2.*pq0(x2);

% partial derivatives
dpq_1 =  b1.*dpq0(x1) + b2.*dpq0(x2);
dpq_2 = -b1.*dpq0(x1) + b2.*dpq0(x2);
dpq_3 = (-pq0(x1) + pq0(x2)).*ones(size(kappa));

% collect all partial derivatives in an array
Dpq = cat(5,dpq_1,dpq_2,dpq_3);

% compute the entries of the CFIM for SPADE
I_HG = zeros([3,3,numel(x0),numel(s),numel(kappa)]);
for i = 1:3
    for j = i:3
        I_HG(i,j,:,:,:) = sum(Dpq(:,:,:,:,i).*Dpq(:,:,:,:,j)./(pq+1e-25),1);
    end
end
I_HG = I_HG + pagetranspose(I_HG.*(1-eye(3)));
I_HG(1:2,1:2,:,:,:) = I_HG(1:2,1:2,:,:,:)*sigma^2; % normalize physical dimensions

% compute the CRB matrix for SPADE
CRB_HG = pageinv(I_HG);

% MATHEMATICA OUTPUT
% Compute the CFIM in the sub-diffraction taylor approximation
I_HG_approx = zeros(size(I_HG));
I_HG_approx(1,1,:,:,:) = (1-(s/sigma).^2.*(1-4*kappa.^2).*(1+(sigma./x0).^2));
I_HG_approx(1,2,:,:,:) = (2*kappa + (s./x0).*(1-4*kappa.^2) - 2*kappa.*(s/sigma).^2 .* (1-4*kappa.^2).*(1+2*(sigma./x0).^2));
I_HG_approx(1,3,:,:,:) = (2*s.*(x0-2*s.*kappa))./x0;
I_HG_approx(2,2,:,:,:) = (4*kappa.^2 + 4*kappa.*s./x0 .*(1-4*kappa.^2) + (s./x0).^2.*(1-4*kappa.^2).*(2*(x0./sigma).^2.*(1-6*kappa.^2)+(1-16*kappa.^2)));
I_HG_approx(2,3,:,:,:) = (2*kappa.*s + 2*s.^2./x0.*(1-4*kappa.^2));
I_HG_approx(3,3,:,:,:) = (4*s.^2).*ones(size(x0)).*ones(size(kappa));
I_HG_approx = I_HG_approx./(sigma^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% FISHER INFO FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% coloring and legend details
kappa_linecolors = hsv(numel(kappa)+2)*.7;
x0_linecolors = hsv(numel(x0)+2)*.7;
x0_legend_names = arrayfun(@(j)sprintf('%.2f',x0(j)/sigma),1:numel(x0),'UniformOutput',false);
kappa_legend_names = arrayfun(@(j)sprintf('%.2f',kappa(j)),1:numel(kappa),'UniformOutput',false);
parameter_names = {'x_0','s','\kappa'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% QUANTUM BOUND %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the diagonal elements of the QFIM
if plot_list(1)
figure
t0 = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for i = 1:3
    nexttile(i)
    hold on
    for k = 1:numel(kappa)
        plot(squeeze(s)/sigma,squeeze(Q(i,i,1,:,k)),'Color',kappa_linecolors(k,:),'LineWidth',1.5)
    end
    hold off
    xlabel('$s/\sigma$','interpreter','latex')
    ylabel(['Fisher Information ','$Q_{',parameter_names{i},'}(s)\cdot\sigma^2$'],'interpreter','latex')
    if i==3
        leg=legend(kappa_legend_names);
        leg.Location='NorthWest';
        title(leg,'$\kappa$','interpreter','latex')
        ylabel(['Fisher Information ','$Q_{',parameter_names{i},'}(s)$'],'interpreter','latex')
        %set(gca,'yscale','log')
    end
    axis square
    box on
    title(['$',parameter_names{i},'$'],'interpreter','latex')
end
title(t0,{'Diagonal Quantum Fisher Information Matrix Elements'},'interpreter','latex')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% DIRECT IMAGING %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_list(2)
% plot the diagonal elements of the CFIM
figure
t1 = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for i = 1:3
    nexttile(i)
    hold on
    for k = 1:numel(kappa)
        plot(squeeze(s)/sigma,squeeze(I_DD(i,i,1,:,k)),'Color',kappa_linecolors(k,:),'LineWidth',1.5)
    end
    hold off
    xlabel('$s/\sigma$','interpreter','latex')
    ylabel(['Fisher Information  ','$I_{',parameter_names{i},'}^{[DI]}(s)\cdot \sigma^2$'],'interpreter','latex')
    if i==3
        leg=legend(kappa_legend_names);
        leg.Location='NorthWest';
        title(leg,'$\kappa$','interpreter','latex')
        ylabel(['Fisher Information  ','$I_{',parameter_names{i},'}^{[DI]}(s)$'],'interpreter','latex')
        %set(gca,'yscale','log')
    end
    axis square
    box on
    title(['$',parameter_names{i},'$'],'interpreter','latex')
end
title(t1,{'Diagonal Fisher Information Matrix Elements','Direct Imaging'},'interpreter','latex')
end

% plot the diagonal elements of the CRB matrix
if plot_list(3)
figure
t2=tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for i = 1:3
    nexttile(i)
    hold on
    for k = 1:numel(kappa)
        plot(squeeze(s)/sigma,squeeze(CRB_DD(i,i,1,:,k)),'Color',kappa_linecolors(k,:),'LineWidth',1.5)
    end
    hold off
    xlabel('$s/\sigma$','interpreter','latex')
    ylabel(['Cramer-Rao Bound ','$\Sigma_{',parameter_names{i},'}^{[DI]}(s)\cdot sigma^2$'],'interpreter','latex')
    if i==3
        leg=legend(kappa_legend_names);
        title(leg,'$\kappa$','interpreter','latex')
        ylabel(['Cramer-Rao Bound ','$\Sigma_{',parameter_names{i},'}^{[DI]}(s)$'],'interpreter','latex')
    end
    set(gca,'yscale','log')
    axis square
    box on
    title(['$',parameter_names{i},'$'],'interpreter','latex')
end
title(t2,{'Diagonal Cramer-Rao Matrix Elements','Direct Imaging'},'interpreter','latex')
end

if plot_list(4)
% Compare approximate (and analytic) CFIM to numerical CFIM 
YLIMS = [.4,1.1; -.1,1; -.1,2;-.1,2;-.1,1.1;-1,10];
figure
t3=tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
n=1;
for i=1:3
    for j=i:3
        %nexttile(sub2ind([3,3],j,i))
        nexttile(n)
        hold on
        for k = 1:numel(kappa)
            plot(squeeze(s)/sigma,squeeze(I_DD(i,j,1,:,k)),'Color',kappa_linecolors(k,:),'LineWidth',1.5);
        end
        for k = 1:numel(kappa)
            plot(squeeze(s)/sigma,squeeze(I_DD_approx(i,j,1,:,k)),'--','Color',kappa_linecolors(k,:),'LineWidth',1.5);
        end
        hold off
        axis square
        box on
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        temp = ylim;
        ylim(YLIMS(n,:))
        xlabel('$s/\sigma$','interpreter','latex')
        title({sprintf('$I^{[DD]}_{%d%d}(s)$',i,j)},'interpreter','latex')
        n=n+1;
    end
    if i==3
        leg=legend(kappa_legend_names);
        title(leg,'$\kappa$','interpreter','latex')
    end
    
end
title(t3,{'Analytic v. Numeric CFIM','Direct Detection'},'interpreter','latex')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% HG-SPADE %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot diagonal elements of the Fisher Information matrix as a function of
% separation for different values of x0 (with kappa=0)
if plot_list(5)
figure
t4 = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for i = 1:3
    nexttile(i)
    hold on
    for k = 1:numel(x0)
        plot(squeeze(s)/sigma,squeeze(I_HG(i,i,k,:,1)),'Color',kappa_linecolors(k,:),'LineWidth',1.5)
    end
    hold off
    xlabel('$s/\sigma$','interpreter','latex')
    ylabel(['Fisher Information  ','$I_{',parameter_names{i},'}^{[HG]}(s)\cdot\sigma^2$'],'interpreter','latex')
    if i==3
        leg=legend(x0_legend_names);
        leg.Location='NorthWest';
        title(leg,'$x_0/\sigma$','interpreter','latex')
        ylabel(['Fisher Information  ','$I_{',parameter_names{i},'}^{[HG]}(s)$'],'interpreter','latex')
        %set(gca,'yscale','log')
    end
    xlim([min(s),max(s)]/sigma)
    axis square
    box on
    title(['$',parameter_names{i},'$'],'interpreter','latex')
end
title(t4,{'Diagonal Fisher Information Matrix Elements','HG SPADE',sprintf('$\\kappa=%.f$',kappa(1))},'interpreter','latex')
end

% plot diagonal elements of the Fisher Information matrix as a function of
% separation for different values of kappa (with x0 = 0)
if plot_list(6)
figure
t5 = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for i = 1:3
    nexttile(i)
    hold on
    for k = 1:numel(kappa)
        plot(squeeze(s)/sigma,squeeze(I_HG(i,i,1,:,k)),'Color',kappa_linecolors(k,:),'LineWidth',1.5)
    end
    hold off
    xlabel('$s/\sigma$','interpreter','latex')
    ylabel(['Fisher Information  ','$I_{',parameter_names{i},'}^{[HG]}(s)/\sigma^2$'],'interpreter','latex')
    if i==3
        leg=legend(kappa_legend_names);
        title(leg,'$\kappa$','interpreter','latex')
        ylabel(['Fisher Information  ','$I_{',parameter_names{i},'}^{[HG]}(s)$'],'interpreter','latex')
        %set(gca,'yscale','log')
    end
    xlim([min(s),max(s)]/sigma)
    ylim([0,1])
    axis square
    box on
    title(['$',parameter_names{i},'$'],'interpreter','latex')
end
title(t5,{'Diagonal Fisher Information Matrix Elements','HG SPADE',sprintf('$x_0 / \\sigma=%.f$',x0(1))},'interpreter','latex')
end

% Plot the diagonal elements of the CRB matrix for different values of
% kappa (with x0 = 0)
if plot_list(7)
figure
t6=tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for i = 1:3
    nexttile(i)
    hold on
    for k = 1:numel(kappa)
        plot(squeeze(s)/sigma,squeeze(CRB_HG(i,i,1,:,k)),'Color',kappa_linecolors(k,:),'LineWidth',1.5)
    end
    hold off
    xlabel('$s/\sigma$','interpreter','latex')
    ylabel(['Cramer-Rao Bound ','$\Sigma_{',parameter_names{i},'}^{[HG]}(s) \sigma^2$'],'interpreter','latex')
    if i==3
        leg=legend(kappa_legend_names);
        title(leg,'$\kappa$','interpreter','latex')
        ylabel(['Cramer-Rao Bound ','$\Sigma_{',parameter_names{i},'}^{[HG]}(s)$'],'interpreter','latex')
    end
    xlim([min(s),max(s)]/sigma)
    set(gca,'yscale','log')
    axis square
    box on
    title(['$',parameter_names{i},'$'],'interpreter','latex')
end
title(t6,{'Diagonal Cramer-Rao Matrix Elements','HG-SPADE',sprintf('$x_0/\\sigma=%.f$',x0(1))},'interpreter','latex')
end

% Plot the diagonal elements of the CRB matrix for different values of x0
% (with kappa = 0)
if plot_list(8)
figure
t7=tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for i = 1:3
    nexttile(i)
    hold on
    for k = 1:numel(x0)
        plot(squeeze(s)/sigma,squeeze(CRB_HG(i,i,k,:,1)),'Color',x0_linecolors(k,:),'LineWidth',1.5)
    end
    hold off
    xlabel('$s/\sigma$','interpreter','latex')
    ylabel(['Cramer-Rao Bound ','$\Sigma_{',parameter_names{i},'}^{[HG]}(s) \sigma^2$'],'interpreter','latex')
    if i==3
        leg=legend(x0_legend_names);
        title(leg,'$x_0/\sigma$','interpreter','latex')
         ylabel(['Cramer-Rao Bound ','$\Sigma_{',parameter_names{i},'}^{[HG]}(s)$'],'interpreter','latex')
    end
    xlim([min(s),max(s)]/sigma)
    set(gca,'yscale','log')
    axis square
    box on
    title(['$',parameter_names{i},'$'],'interpreter','latex')
end
title(t7,{'Diagonal Cramer-Rao Matrix Elements','HG-SPADE',sprintf('$\\kappa=%.f$',kappa(1))},'interpreter','latex')
end

% Compare approximate (and analytic) CFIM to numerical CFIM for HG SPADE
if plot_list(9)
figure
t8=tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
n=1;
for i=1:3
    for j=i:3
        %nexttile(sub2ind([3,3],j,i))
        nexttile(n)
        hold on
        for k = 1:numel(x0)
            plot(squeeze(s)/sigma,squeeze(I_HG(i,j,k,:,1)),'Color',x0_linecolors(k,:),'LineWidth',1.5);
        end
        for k = 1:numel(x0)
            plot(squeeze(s)/sigma,squeeze(I_HG_approx(i,j,k,:,1)),'--','Color',x0_linecolors(k,:),'LineWidth',1.5);
        end
        hold off
        axis square
        box on
        set(gca,'xscale','log')
        temp = ylim;
        ylim(YLIMS(n,:))
        xlabel('$s/\sigma$','interpreter','latex')
        title({sprintf('$I^{[HG]}_{%d%d}(s)$',i,j)},'interpreter','latex')
        set(gca,'yscale','log')
        n=n+1;
    end
    if i==3
        leg=legend(x0_legend_names);
        title(leg,'$x_0/\sigma$','interpreter','latex')
       
    end
    
end
title(t8,{'Analytic v. Numeric CFIM','HG SPADE',sprintf('$\\kappa=%f.$',kappa(1))},'interpreter','latex')
end