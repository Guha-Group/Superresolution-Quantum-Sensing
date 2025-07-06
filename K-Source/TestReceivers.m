% A script for testing receivers dedicated to multi-source color center
% sensing.
% 
% Author: Nico Deshler
% Date: 4/28/2025

% add utils path
addpath('utils/')

% ground truth parameters
sigma = 1;                          
num_sources = 5;
d_min = sigma/4;
contrast = .1;
visualize_flag = 0;

% generate scene
xyb = GenerateRandomConstellation(num_sources,d_min,contrast,sigma);
%xy = min_sep_frac*rand([num_sources,2]); xy = xy-mean(xy);
%b = dirchrnd(ones(num_sources,1)); xyb = [xy,b];
[~,srt_id] = sort(xyb(:,3),1); xyb = xyb(srt_id,:);

% Photon Allocations
M = 1e6;           % number of calibration photons
N = 1e6;           % number of sensing photons

% simulate SPADE estimation receiver
splitting_ratio = .1;
n_max = 10;
[n,m] = HGIndices(n_max);
[xyb_SS, out_SS] = SimulateReceiver(xyb,M,N,'StaticSPADE',splitting_ratio,[n,m],sigma,visualize_flag);

% simulate direct imaging receiver
[xyb_DI, out_DI] = SimulateReceiver(xyb,M,N,'DirectImaging',sigma,visualize_flag);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT RECONSTRUCTIONS FROM EACH RECEIVER

% reference circle for diffraction limit
theta = linspace(0,2*pi,1001);
radius = .5*ones(size(theta));
[xc,yc] = pol2cart(theta,radius);


figure;
T = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
title(T,'Receiver Performance Comparison','interpreter','latex')

spc = [0.2,0.8,1];
nexttile(1)
hold on
scatter(xyb(:,1)/sigma,xyb(:,2)/sigma,'k','filled')
scatter(xyb_SS(:,1)/sigma,xyb_SS(:,2)/sigma,'d','MarkerFaceColor',spc,'MarkerFaceAlpha',0,'MarkerEdgeColor',spc,'LineWidth',1.5)
scatter(xyb_DI(:,1)/sigma,xyb_DI(:,2)/sigma,'r','s','MarkerFaceAlpha',.5,'MarkerEdgeColor','r','LineWidth',1.5)
plot(xc,yc,'--k','LineWidth',.5)
scatter(xyb_SS(:,1)/sigma,xyb_SS(:,2)/sigma,'d','MarkerFaceColor',spc,'MarkerFaceAlpha',0,'MarkerEdgeColor',spc,'LineWidth',1.5)
hold off
xlim(.5*[-1,1]); ylim(.5*[-1,1])
xticks(-1:.25:1); yticks(-1:.25:1);
axis square; box on;
legend({'Ground Truth','SPADE','Direct Imaging','Diffraction Spot',''},'interpreter','latex')
xlabel('$x/\sigma$','interpreter','latex'); ylabel('$y/\sigma$','interpreter','latex');
title('Position Estimates','interpreter','latex');

nexttile(2)
barArray = bar((1:num_sources).',[xyb(:,3),sort(xyb_SS(:,3)),sort(xyb_DI(:,3))],'FaceColor','flat');
barArray(1).CData = [0,0,0]; barArray(1).FaceAlpha=.5;
barArray(2).CData = spc; barArray(2).FaceAlpha=.5;
barArray(3).CData = [1,0,0]; barArray(3).FaceAlpha=.5;
axis square; box on;
title('Brightness Estimates','interpreter','latex')
xlabel('Source Index $k$','interpreter','latex')
ylabel('Relative Brightness $b_k$','interpreter','latex')
legend({'Ground Truth','SPADE','Direct Imaging'},'interpreter','latex')


% plot the YKL modes
[X,Y] = meshgrid(sigma*linspace(-5,5,1001));
PSF_modes = 1/sqrt(2*pi*sigma^2)*...
            exp(- (1/(2*sigma)^2)*(...
            (X-permute(xyb_SS(:,1),[3,2,1])).^2 +...
            (Y-permute(xyb_SS(:,2),[3,2,1])).^2 ));

figure
T = tiledlayout(floor(sqrt(num_sources)),ceil(sqrt(num_sources)),'TileSpacing','compact','Padding','compact');
title(T,'YKL/Helstrom Modes','interpreter','latex')
A = out_SS.YKL_GaussRepn;
for t=1:num_sources
    nexttile(t)
    mode = sum(permute(A(t,:),[3,1,2]).*PSF_modes,3);
    imagescComplex([min(X(:)),max(X(:))]/sigma,[min(Y(:)),max(Y(:))]/sigma,mode,.7,0);
    %colormap(turbo)
    %imagesc([min(X(:)),max(X(:))]/sigma,[min(Y(:)),max(Y(:))]/sigma,abs(mode).^2);
    axis square
    axis off
    set(gca,'ydir','normal')
    xlabel('$x/\sigma$','interpreter','latex')
    ylabel('$y/\sigma$','interpreter','latex')
    title(sprintf('$|\\pi_{%d}\\rangle$',t),'interpreter','latex')

end





%{

nexttile(3)
hold on
scatter(xyb(:,1)/sigma,xyb(:,2)/sigma,30*num_sources*xyb(:,3),'k','filled')
scatter(xyb_SS(:,1)/sigma,xyb_SS(:,2)/sigma,30*num_sources*xyb_SS(:,3),'b','d','filled','MarkerFaceAlpha',.5,'MarkerEdgeColor','b','LineWidth',1.5)
scatter(xyb_DI(:,1)/sigma,xyb_DI(:,2)/sigma,30*num_sources*xyb_DI(:,3),'r','s','filled','MarkerFaceAlpha',.5,'MarkerEdgeColor','r','LineWidth',1.5)
plot(xc,yc,'--k','LineWidth',.5)
hold off
xlim([-1,1]); ylim([-1,1])
xticks(-1:.5:1); yticks(-1:.25:1);
axis square; box on;
legend({'Ground Truth','SPADE','Direct Imaging','Diffraction Spot'},'interpreter','latex')
xlabel('$x/\sigma$','interpreter','latex'); ylabel('$y/\sigma$','interpreter','latex');
title('Reconstruction','interpreter','latex');

%}

% Auxiliary figure with numbering
figure
hold on
scatter(xyb(:,1)/sigma,xyb(:,2)/sigma,'k','filled')
scatter(xyb_SS(:,1)/sigma,xyb_SS(:,2)/sigma,'d','MarkerFaceColor',spc,'MarkerFaceAlpha',0,'MarkerEdgeColor',spc,'LineWidth',1.5)
scatter(xyb_DI(:,1)/sigma,xyb_DI(:,2)/sigma,'r','s','MarkerFaceAlpha',.5,'MarkerEdgeColor','r','LineWidth',1.5)
plot(xc,yc,'--k','LineWidth',.5)
scatter(xyb_SS(:,1)/sigma,xyb_SS(:,2)/sigma,'d','MarkerFaceColor',spc,'MarkerFaceAlpha',0,'MarkerEdgeColor',spc,'LineWidth',1.5)
hold off
xlim(.5*[-1,1]); ylim(.5*[-1,1])
xticks(-1:.25:1); yticks(-1:.25:1);
axis square; box on;
legend({'Ground Truth','SPADE','Direct Imaging','Diffraction Spot',''},'interpreter','latex')
xlabel('$x/\sigma$','interpreter','latex'); ylabel('$y/\sigma$','interpreter','latex');
title('Position Estimates','interpreter','latex');

hold on
for k = 1:num_sources
    text(xyb(k,1)/sigma+.01,xyb(k,2)/sigma+.01, 0, num2str(k))
end
hold off
