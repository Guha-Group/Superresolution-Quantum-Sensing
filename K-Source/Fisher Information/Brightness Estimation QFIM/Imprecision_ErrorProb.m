% constellation params
num_sources = 3;
num_constellations = 5;

% make the constellations
r = ones(num_sources,1);
th = pi/2 + 2*pi/3 * (0:(num_sources-1))';
[x,y] = pol2cart(th,r);
disp_x = zeros(num_sources,1,num_constellations);
disp_x(1,1,:) = (0:(num_constellations-1))/2;
disp_y = zeros(size(disp_x));
xx = x + disp_x;
yy = y + disp_y;
xy_constellations = [xx,yy];
%{
figure
hold on
for k = 1:5
    fill(xy_constellations(:,1,k),xy_constellations(:,2,k),'b','FaceAlpha',.2);
end
hold off
%}
%xy_constellations = 2*(rand(num_sources,2,num_constellations)-.5);
%xy_constellations = xy_constellations - mean(xy_constellations,1);

% generate simplex samples
num_samples = 50;
[simplex_samples,num_simplex_samples] = utils.sample_3_simplex(num_samples);

% precision and error probability containers
H = zeros(num_simplex_samples,num_constellations);
P = zeros(num_simplex_samples,num_constellations);
for n = 1:num_constellations
    
    % get the constellation
    xy = xy_constellations(:,:,n);
    
    % loop through simplex samples computing the QFIM for param estimation
    % and the P_e_min for state discrimination
    QFIM_stack = zeros(num_sources-1,num_sources-1,num_simplex_samples,num_sources);
    P_e_min = zeros(num_simplex_samples,1);
    for k = 1:num_simplex_samples
        
        % get the states and priors associated with the scene
        xyb = [xy,simplex_samples(:,k)];
        states = utils.scene_states(xyb);
        priors = xyb(:,3);
    
        % compute QFI matrix for the brightness parameters
        for t = 1:num_sources % run through all indexing permutations
            states = circshift(states,1,2);
            priors = circshift(priors,1,1);
            QFIM_stack(:,:,k,t) = utils.brightness_params_QFIM(states, priors);
        end
    
        % compute the min probability of error for state discrimination
        [~,p_e] = utils.ykl_states(states, priors);
        P_e_min(k) = p_e;
        sprintf('(%d,%d) of (%d,%d)',n,k,num_constellations,num_simplex_samples)
    end
    
    % compute the estimation imprecision and the minimum probability of
    % error for state discrimination
    H(:,n) = squeeze(mean(sum(eye(num_sources-1).*pageinv(QFIM_stack),[1,2]),4)); % imprecision
    P(:,n) = P_e_min;
end
%% FIGURES
%}
% normalize the metrics
%h = H./max(H(:));
%p = P./max(P(:));
%h = H./max(H,[],1);
%p = P./max(P,[],1);
h = H;
p = P;

% reference circle for Rayleigh limit
[xr,yr]=pol2cart(linspace(0,2*pi,100),ones(1,100));

figure
N = num_constellations-1;
tiledlayout(2,(N),'TileSpacing','compact','Padding','compact')

for n = 1:N

    x = xy_constellations(:,1,n);
    y = xy_constellations(:,2,n);
    
    nexttile(n)
    plot_simplex_heatmap(x,y,simplex_samples,h(:,n))
    clim([0,max(h(:))])
    hold on
    axis square
    box on
    scatter(x,y,10,'k','filled')
    xlim(1.5*[-1,1]); ylim(1.5*[-1,1]);
    xlabel('$x/\sigma$','interpreter','latex')
    ylabel('$y/\sigma$','interpreter','latex')
    %title({'Brightness Estimation','Imprecision'},'interpreter','latex')
    if n==N
        cbar=colorbar;
        ylabel(cbar,'Min. Imprecision','interpreter','latex');
    end
    hold off

    % plot the interior point for the probability of error
    nexttile(n+N)
    plot_simplex_heatmap(x,y,simplex_samples,p(:,n))
    clim([0,max(p(:))])
    hold on
    axis square
    box on
    scatter(x,y,10,'k','filled')
    xlim(1.5*[-1,1]); ylim(1.5*[-1,1]);
    xlabel('$x/\sigma$','interpreter','latex')
    ylabel('$y/\sigma$','interpreter','latex')
    %title({'State Classification','Error Probability'},'interpreter','latex')
    if n==N
        cbar=colorbar;
        ylabel(cbar,'Min. Error Probability','interpreter','latex');
    end
    hold off
end


function plot_simplex_heatmap(x, y, simplex_samples, value)
% plots a heatmap of the 'value' of the simplex samples defined by the
% vertices at x,y
    
    % compute interior points of the simplex with barycentric coordinates
    xb = simplex_samples'*x;
    yb = simplex_samples'*y;
    
    % reference circle for Rayleigh limit
    [xr,yr]=pol2cart(linspace(0,2*pi,100),ones(1,100));

    % meshgrid for the simplex support
    [X,Y] = meshgrid(linspace(min([x;y]),max([x;y]),500));
    
    % interpolate the function over the simplex
    F = scatteredInterpolant(xb,yb,value,'linear');
    Z = F(X,Y);
    
    % remove interpolated points that lie outside the simplex
    inpts = inpolygon(X(:),Y(:),x,y);
    Z(~inpts) = nan;

    % plot the function over the simplex
    surf(X,Y,Z)
    colormap(turbo)
    shading interp
    view(0,90)
    hold on
    simplex = polyshape(x,y);
    plot(simplex,'FaceAlpha',0,'LineWidth',4)
    plot(xr,yr,'--k','LineWidth',1)
    hold off
    axis square
    grid on
    xlim(1.5*[-1,1])
    ylim(1.5*[-1,1])
end
