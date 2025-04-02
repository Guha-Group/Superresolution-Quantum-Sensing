sigma = 1;
num_sources = 3;

% generate constellation
r = 1*ones(1,3);
theta = 2*pi*rand(1,3);
[sx,sy] = pol2cart(theta,r);
%sx = sigma*(3*(rand(1,3)-.5));
%sy = sigma*(3*(rand(1,3)-.5));
sx = sx-mean(sx); sy = sy-mean(sy);
sources_xy = sigma*[sx',sy'];

% generate all combinations of possible weight vectors on the tri simplex via
% 'stick-breaking'
q = linspace(0,1,50);
num_simplex_samples = (numel(q)+1)*numel(q)/2;
simplex_samples = zeros(3,num_simplex_samples);
k=1;
for i=1:numel(q)
    for j=i:numel(q)
        p1 = q(i); p2 = q(j)-q(i); p3 = 1-p1-p2;
        simplex_samples(:,k) = [p1,p2,p3]';
        k = k+1;
    end
end
% prune samples that lie on boundary of simplex
simplex_samples = simplex_samples(:,prod(simplex_samples,1)>0);
num_simplex_samples = size(simplex_samples,2);

% evaluate the trace of the  QFIM matrix for all sample points on the simplex
QFIM_stack = zeros([num_sources,num_sources,num_simplex_samples]);
for k=1:size(simplex_samples,2)
    QFIM_stack(:,:,k) = QFIM(simplex_samples(:,k),sources_xy);
end
QCRB_stack = pageinv(QFIM_stack);
QCRBtr = squeeze(sum(QCRB_stack.*eye(num_sources),[1,2]));

% precision
H = 1./QCRBtr;
% remove outliers
H(H<0) = nan;
H(H>5) = nan;
H = H/max(H);


% plot the precision over the simplex
figure
sx = sources_xy(:,1)'; sy= sources_xy(:,2)';
simplex = polyshape(sx,sy);
plot(simplex,'FaceColor','k','FaceAlpha',0.1)
hold on
colors = winter(num_simplex_samples+1);
for k=1:num_simplex_samples
    if ~isnan(H(k))
    scatter(sx*simplex_samples(:,k),...
        sy*simplex_samples(:,k),'filled',...
        'MarkerFaceColor',colors(floor(H(k)*num_simplex_samples +1),:)*.9,...
        'MarkerFaceAlpha',H(k)*.9+.1)
    end
end
hold off
axis equal

xlim([-1,1])
ylim([-1,1])





% max normed results for plotting colors
QCRBtr_maxnorm = abs(QCRBtr);
QCRBtr_maxnorm(QCRBtr_maxnorm<0)=nan;
QCRBtr_maxnorm = log(abs(QCRBtr_maxnorm));
QCRBtr_maxnorm = QCRBtr_maxnorm - min(QCRBtr_maxnorm);
QCRBtr_maxnorm = QCRBtr_maxnorm/max(QCRBtr_maxnorm);

% parametric curve for point source locations



% plot the estimates over the simplex relative to the ground truth
figure
sx=sources_xy(:,1)'; sy= sources_xy(:,2)';
simplex = polyshape(sx,sy);
plot(simplex,'FaceColor','k','FaceAlpha',0.1)
hold on
colors = winter(num_simplex_samples+1);
for k=1:num_simplex_samples
    if ~isnan(QCRBtr_maxnorm(k))
    scatter(sx*simplex_samples(:,k),...
        sy*simplex_samples(:,k),'filled',...
        'MarkerFaceColor',colors(floor((1-QCRBtr_maxnorm(k))*num_simplex_samples +1),:)*.9,...
        'MarkerFaceAlpha',1-QCRBtr_maxnorm(k))
    end
end
hold off
axis equal



% plot the simplex as a surface
figure
F = scatteredInterpolant((sx*simplex_samples)',(sy*simplex_samples)',QCRBtr_maxnorm,'linear');
[X,Y] = meshgrid(linspace(min([sx,sy]),max([sx,sy]),200));
Z = F(X,Y);
inpts = inpolygon(X(:),Y(:),sx,sy);
Z(~inpts) = nan;
surf(X,Y,-Z)
colormap(winter)
view(0,90)
shading interp
hold on
plot(simplex,'FaceAlpha',0,'LineWidth',4)
[xc,yc]=pol2cart(linspace(0,2*pi,100),ones(1,100));
plot(xc,yc,'--k','LineWidth',1.5)
hold off
xlim(1*[-1,1])
ylim(1*[-1,1])
axis square
cbar = colorbar;
cbar.Ticks = [min(cbar.Ticks),max(cbar.Ticks)];
cbar.TickLabels = {'Imprecise','Precise'};
cbar.TickLabelInterpreter = 'latex';
ylabel(cbar,'-LOG(Tr$[\Sigma_{Q}(\mathbf{p})])$  Arb. Units','interpreter','latex')
xlabel('$x/\sigma$','interpreter','latex')
ylabel('$y/\sigma$','interpreter','latex')
yticks(xticks)

function Q = QFIM(weights, sources_xy)
% computes the QFIM matrix for the compositional parameters of a
% constellation with positions sources_xy;

% get eigenvalues and their numeric derivatives
[R,L] =  EigenvectorRepresentationMatrix(weights,sources_xy);
[dR,dL] = DerivativeEigRepMat(weights,sources_xy);

% get inner product matrices for derivative eigenstates
G = (inv(R))'*inv(R);
N  = (pagemtimes(pagemtimes(pagectranspose(dR),G),dR).^(-1/2)).*eye(size(G,1));
B = pagemtimes(inv(R),pagemtimes(dR,N));


% compute the QFIM
num_sources = numel(weights);
Q=zeros(num_sources);
for i=1:num_sources
    for j =i:num_sources
        T1 = dL(:,i).*dL(:,j)./L;
        T1 = sum(T1(L>0));
        T2 = 2*(L-L.').^2./(L+L.').*abs(real(B(:,:,i).*(B(:,:,j)')));
        T2 = sum(T2((L+L.')>0),'all');
        Q(i,j) = T1 + T2;
    end
end
Q = Q+(Q.'.*(1-eye(num_sources)));
end



function [R,L] = EigenvectorRepresentationMatrix(weights,sources_xy)
    % R is a matrix that maps from the non-orthogonal (yet linearly
    % independent) states |psi> to the orthogonal eigenvectors of the
    % density operator |lambda>. That is |lambda> = R|psi>.
    %
    % L is a vector containing the eigenvalues of the density
    % operator.
    
    % construct some helpful structures
    P = diag(weights);
    deltas_xy = permute(sources_xy,[1,3,2])-permute(sources_xy,[3,1,2]);
    
    % Gram Matrix
    G = exp(-sum(deltas_xy.^2,3)/8);
    
    % spectral decomposition of the Gram Matrix
    [U,D] = eig(G);
    
    % Evaluate S matrix
    A = sqrt(P)*U*sqrt(D);
    S = A'*A;
    [V_dagger,L] = eig(S,'vector');
    V = V_dagger';
    
    % Representation matrix
    R = U*diag(diag(D).^(-1/2))*V';
end

function [dR,dL] = DerivativeEigRepMat(weights,sources_xy)
    num_sources = numel(weights);
    dR = zeros(num_sources,num_sources,num_sources);
    dL = zeros(num_sources,num_sources);
    epsilon = 1e-14;
    for k = 1:num_sources
        [R1,L1] = EigenvectorRepresentationMatrix(weights,sources_xy);
        %dp = zeros(num_sources,1); dp(k) = epsilon;
        dp = - epsilon*ones(num_sources,1); dp(k) = dp(k) + epsilon;
        weights_dp = weights+dp; 
        weights_dp(weights_dp<0) = 0;
        weights_dp(weights_dp>1) = 1;
        weights_dp = weights_dp/sum(weights_dp);
        [R2,L2] = EigenvectorRepresentationMatrix(weights_dp,sources_xy);
        dR(:,:,k) = (R2-R1)./epsilon;
        dL(:,k) = (L2-L1)./epsilon; 
    end

end