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
d_min = sigma/20;
contrast = .1;
visualize_flag = 0;

% generate scene
xyb = GenerateRandomConstellation(num_sources,d_min,contrast,sigma);
%xy = min_sep_frac*rand([num_sources,2]); xy = xy-mean(xy);
%b = dirchrnd(ones(num_sources,1)); xyb = [xy,b];
[~,srt_id] = sort(xyb(:,3),1); xyb = xyb(srt_id,:);

% Photon Allocations
M = 1e7;           % number of calibration photons
N = 1e7;           % number of sensing photons

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

receiver_colors = .9*[75,212,219;      % Direct Imaging
                   72,219,82]/256;  % SPADE


figure;
T = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
title(T,'Receiver Performance Comparison','interpreter','latex')

nexttile(1)
hold on
scatter(xyb(:,1)/sigma,xyb(:,2)/sigma,'k','filled')
scatter(xyb_DI(:,1)/sigma,xyb_DI(:,2)/sigma,'s','MarkerFaceColor',receiver_colors(1,:),'MarkerFaceAlpha',.5,'MarkerEdgeColor',receiver_colors(1,:),'LineWidth',1.5)
scatter(xyb_SS(:,1)/sigma,xyb_SS(:,2)/sigma,'d','MarkerFaceColor',receiver_colors(2,:),'MarkerFaceAlpha',.5,'MarkerEdgeColor',receiver_colors(2,:),'LineWidth',1.5)
plot(xc,yc,'--k','LineWidth',.5)
hold off
xlim(.5*[-1,1]); ylim(.5*[-1,1])
xticks(-1:.25:1); yticks(-1:.25:1);
axis square; box on;
legend({'Ground Truth','Direct Imaging','SPADE','Diffraction Spot'},'interpreter','latex','Location','Northwest')
xlabel('$x/\sigma$','interpreter','latex'); ylabel('$y/\sigma$','interpreter','latex');
title('Position Estimates','interpreter','latex');
hold on
for k = 1:num_sources
    text(xyb(k,1)/sigma+.02,xyb(k,2)/sigma+.02, 0, num2str(k))
end
text(.3*(sigma/2),.8*(sigma/2), 0,sprintf('$d_{\\mathrm{min}} = %.2f \\sigma$',d_min/sigma),'interpreter','latex')
hold off


nexttile(2)
%barArray = bar((1:num_sources).',[xyb(:,3),sort(xyb_DI(:,3)),sort(xyb_SS(:,3))],'FaceColor','flat');
barArray = bar((1:num_sources).',[xyb(:,3),xyb_DI(:,3),xyb_SS(:,3)],'FaceColor','flat');
barArray(1).CData = [0,0,0]; barArray(1).FaceAlpha=.5;
barArray(2).CData = receiver_colors(1,:); barArray(2).FaceAlpha=.5;
barArray(3).CData = receiver_colors(2,:); barArray(3).FaceAlpha=.5;
axis square; box on;
title('Brightness Estimates','interpreter','latex')
xlabel('Source Index $k$','interpreter','latex')
ylabel('Relative Brightness $b_k$','interpreter','latex')
legend({'Ground Truth','Direct Imaging','SPADE'},'interpreter','latex')

% plot the YKL modes
[X,Y] = meshgrid(sigma*linspace(-5,5,1001));
psfm = 1/sqrt(2*pi*sigma^2)*...
            exp(- (1/(2*sigma)^2)*(...
            (X-permute(xyb_SS(:,1),[3,2,1])).^2 +...
            (Y-permute(xyb_SS(:,2),[3,2,1])).^2 ));


% plot the ykl modes in postion representation
A = out_SS.YKL_GaussRepn;
A = A.*conj(exp(1i*angle(A(:,1))));
ykl_xy = sum(permute(A,[3,4,1,2]).*permute(psfm,[1,2,4,3]),4);

figure
%T = tiledlayout(floor(sqrt(num_sources)),ceil(sqrt(num_sources)),'TileSpacing','compact','Padding','compact');
T = tiledlayout(num_sources,1,'TileSpacing','compact','Padding','compact');
title(T,'YKL Modes','interpreter','latex')
for k=1:num_sources
    nexttile(k)
    xlimits=[min(X(:)),max(X(:))];
    ylimits=[min(Y(:)),max(Y(:))];
    hold on
    imagescComplex(xlimits,ylimits,ykl_xy(:,:,k),.7,0)
    %scatter(xyb_SS(:,1),xyb_SS(:,2),'w','filled')
    plot(xc,yc,'w')
    hold off
    axis square; set(gca,'ydir','normal');
    axis off;
    xlim(xlimits); ylim(ylimits);
    xlabel('$x/\sigma$','interpreter','latex')
    ylabel('$y/\sigma$','interpreter','latex')
    %title(sprintf('$|\\upsilon_{%d}\\rangle$',k),'interpreter','latex')
end