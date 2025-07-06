% Simulates an instance of estimating field parameters from a 
% a modeled brightness response for sub-diffraction color centers.
%
% Available Models:
% -Rabi Oscillation frequency estimation
% -ODMR detuning estimation

addpath('utils/')

% simulation parameters
visualize_flag = 0;         % toggle on/off to see outputs
sigma = 1;                  % Rayleigh unit for gaussian PSF
num_sources = 3;            % number of sources
d_min = sigma/8;            % minimum pair-wise source separation

% pick a model
model_name = 'Rabi';    % {'Rabi','ODMR'};

% HG mode sorter parameters
n_max = 10;                 % max HG index
[n,m] = HGIndices(n_max);   % HG indices
nm = [n,m];                 

% splitting fraction between DI and SPADE for Calibration and Sensing stages
splitting_ratio = .1;       

% generate the locations of the emitters
temp = GenerateRandomConstellationAlt(num_sources,d_min,ones(num_sources,1),sigma);
xy = temp(:,1:2);
xy = xy - mean(xy,1);

% setup brightness variation for the model
switch model_name
    case 'Rabi'
        t = linspace(0,2*pi,200);             % observation time points
        omega_0 = 1;                          % Raw Rabi freq.    
        omega_k = [1,sqrt(2),sqrt(3)]';                   % Rabi frequencies
        contrast = .5;
        %------------------%
        model = @(a,x) RabiModel(x,a,omega_0,contrast);
        params = omega_k;
        param_range = [0,4];

    case 'ODMR'
        t = linspace(0,6,200);               % microwave scan [in linewidths]
        omega_0 = 0;                          % zero-field splitting frequency [in linewidths]
        omega_k = [2,3,4]'; % ODMR detuning frequencies [in linewidths];
        contrast = .5;
        %------------------%
        model = @(a,x) ODMRModel(x,a,omega_0,contrast);
        params = omega_k;
        param_range = [0,5];
end

% get brightness modulation of each source
It = model(params,t);

% relative brightness
bt = It./sum(It,1);

% calibration brightness
b0 = ones(num_sources,1)/num_sources;

% photon allocation
eta0 = 1e6; % emission rate per sample period per source at equal brightness

%M = 1e6; M1 = ceil(splitting_ratio*M); M2 = M-M1;
%N = 1e6; N1 = ceil(splitting_ratio*N); N2 = N-N1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Calibration Stage %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%    
%%% DIRECT IMAGING %%%%
% Localization with direct imaging
M = poissrnd(eta0*num_sources);
xy_samples_DI = DirectImagingMeasurement([xy,b0],sigma,M);
xy_est_DI = DILocalizeSources(xy,xy_samples_DI,sigma);
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SPADE %%%%%%%
% Localization with SPADE
M1 = ceil(splitting_ratio*M);
xy_samples_SS =  DirectImagingMeasurement([xy,b0],sigma,M1);

% centroid estimate
centroid_est = mean(xy_samples_SS,1);

% change pointing 
xy_centered = xy - centroid_est;
xyb_SS = [xy_centered,b0];

% SPADE measurement
M2 = M-M1;
rot_angles = [0,pi/4];
nm_samples = HGSPADEMeasurement(xyb_SS,nm,M2,sigma,rot_angles);

% estimate the source positions from the spade measurement
xy_est_SS = HGLocalizeSources(xy_centered,nm,nm_samples,sigma,rot_angles,visualize_flag);
xy_est_SS = xy_est_SS + centroid_est;
%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
scatter(xy(:,1),xy(:,2),'k','filled')
scatter(xy_est_DI(:,1),xy_est_DI(:,2),'r','filled')
scatter(xy_est_SS(:,1),xy_est_SS(:,2),'b','filled')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Sensing Stage %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% containers for brightness estimates
Nt = zeros(size(t));
b_est_DI = zeros(size(bt));
b_est_SS = zeros(size(bt));
for k = 1:numel(t)

    %%%%%%%%%%%%%%%%%%%%%%%%    
    %%% DIRECT IMAGING %%%%
    % Brightness estimation with Direct Imaging
    N = poissrnd(eta0*sum(It(:,k))); 
    Nt(k) = N; 
    xy_samples_DI = DirectImagingMeasurement([xy,bt(:,k)],sigma,N);
    b_est_DI(:,k) = EstimateBrightnessDI(xy_samples_DI,xy_est_DI,sigma); 
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% SPADE %%%%%%%%
    % Brightness estimation with SPADE
    % direct imaging measurement to get brightness pre-estimates
    N1 = ceil(splitting_ratio*N);
    xy_samples_SS = DirectImagingMeasurement([xy,bt(:,k)],sigma,N1);
    b_pre_est = EstimateBrightnessDI(xy_samples_SS,xy_est_SS,sigma);

    % formulate preliminary scene estimate
    xyb_pre_est = [xy_est_SS,b_pre_est];
    %xyb_pre_est = [xy_est_SS,b0];
    
    % formulate YKL measurement
    N2 = N-N1;
    [YKL_samples,YKL,Psi_est,YKL_GaussRepn] = YKLMeasurement([xy,bt(:,k)],xyb_pre_est,sigma,N2,visualize_flag);

    % estimate the source brightnesses
    b_est_SS(:,k) = EstimateBrightnessYKL(YKL_samples,YKL,Psi_est);
    %%%%%%%%%%%%%%%%%%%%%%%%
end

% use the brightness estimates to recover the target field
It_est_DI = Nt/(eta0).*b_est_DI;
It_est_SS = Nt/(eta0).*b_est_SS;

% shift reconstructed intensities to match models and assign xlabels
switch model_name
    case 'Rabi'
        It_est_DI = It_est_DI - min(It_est_DI,[],2);
        It_est_SS = It_est_SS-min(It_est_SS,[],2);
        xlabel_text = '$t$';
        title_text = 'Rabi Oscillations';
        yrange = [0,1];
    case 'ODMR'
        It_est_SS = It_est_SS-(mean(It_est_SS(:,1),2)-It(:,1));
        xlabel_text = '$(\omega - \omega_0)/\Omega$';
        title_text = 'ODMR Spectrum';
        yrange = [0.7,1];
end

% fit the measurements to the model parameter
param_est_DI = zeros(num_sources,1);
param_est_SS = zeros(num_sources,1);
for j = 1:num_sources
    param_est_DI(j) = ModelFit(t,It_est_DI(j,:),model,param_range);
    param_est_SS(j) = ModelFit(t,It_est_SS(j,:),model,param_range);
end


%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% FIGURES %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
figure;
%-------------------------------------
% Plot the reconstructed constellation
%-------------------------------------
% reference circle for diffraction limit
theta = linspace(0,2*pi,1001);
radius = .5*ones(size(theta));
[xc,yc] = pol2cart(theta,radius);
spc = [0.2,0.8,1];

% plot constellation
hold on
scatter(xy(:,1),xy(:,2),30,'o','k','filled')
scatter(xy_est_SS(:,1),xy_est_SS(:,2),30,'d','MarkerFaceColor',spc,'MarkerFaceAlpha',0,'MarkerEdgeColor',spc,'LineWidth',1.5)
scatter(xy_est_DI(:,1),xy_est_DI(:,2),30,'r','s','MarkerFaceAlpha',.5,'MarkerEdgeColor','r','LineWidth',1.5)
plot(xc,yc,'--k','LineWidth',.5)
hold off
xlim(.5*[-1,1]); ylim(.5*[-1,1])
xticks(-1:.25:1); yticks(-1:.25:1);
axis square; box on;
xlabel('$x/\sigma$','interpreter','latex')
ylabel('$y/\sigma$','interpreter','latex')
title('Source Localization','interpreter','latex')
legend({'Ground Truth','PAD-SPADE','Direct Imaging','Diffraction Spot'},'interpreter','latex')

figure
T = tiledlayout(2,num_sources,'TileSpacing','compact','Padding','compact');
title(T,title_text,'interpreter','latex')


%--------------------------------------
% Plot the brightness estimates and model fit parameters
%--------------------------------------
for k = 1:num_sources
    % Plot the absolute brightness estimates for SPADE
    nexttile(k)
    hold on
    plot(t,It(k,:),'k','LineWidth',1.5)
    plot(t,model(param_est_SS(k),t),'--','Color',spc,'LineWidth',1.5)
    scatter(t,It_est_SS(k,:),20,'filled','MarkerFaceColor',spc,'MarkerFaceAlpha',.3)
    hold off
    axis square; box on;
    ylim(yrange);
    title(sprintf('Source %d',k),'interpreter','latex')
    if k==1
        ylabel({'Source Intensity',sprintf('$I_{%d}(t)$',k)},'interpreter','latex')
    else
        ylabel(sprintf('$I_{%d}(t)$',k),'interpreter','latex')
    end
    if k==num_sources
        legend({'Ground Truth','YKL-SPADE \n (model fit)','YKL-SPADE\n (estimates)'},'interpreter','latex')
    end

    % Plot the absolute brightness estimates for Direct Imaging
    nexttile(k+num_sources)
    hold on
    plot(t,It(k,:),'k','LineWidth',1.5)
    plot(t,model(param_est_DI(k),t),'Color','r','LineWidth',1.5)
    scatter(t,It_est_DI(k,:),20,'filled','MarkerFaceColor','r','MarkerFaceAlpha',.3)
    hold off
    axis square; box on;
    ylim(yrange);
    xlabel(xlabel_text,'interpreter','latex')
    if k==1
        ylabel({'Source Intensity',sprintf('$I_{%d}(t)$',k)},'interpreter','latex')
    else
        ylabel(sprintf('$I_{%d}(t)$',k),'interpreter','latex')
    end
    if k==num_sources
        legend({'Ground Truth','DI (model fit)','DI (estimates)'},'interpreter','latex')
    end
end

%% Functions
function It = RabiModel(t,omega_k,omega_0,contrast)
% returns the optical response of an emitter under Rabi oscillation
It = (omega_0./omega_k).^2 .* (cos(omega_k.*t)).^2; % Rabi emission rate
end

function It = ODMRModel(omega,omega_k,omega_0,contrast)
% Returns the optical response of an emitter under ODMR modulation

% Define Lorentzian and its derivative wrt to the detuning
Lorentzian = @(omega, omega_d) 1 ./ (1 + (omega - omega_0 - omega_d).^2);

% flux rates from either NV source
It = (1 - contrast/2 * (Lorentzian(omega,-omega_k) + Lorentzian(omega,+omega_k)));
end

function a_est = ModelFit(x,y,model,a_range)
    % fits a single parameter model to the data using an iterative refined
    % search for least squares parameter fit
    %
    % assume model is a function handle of the form y = model(a,x) where x
    % is the independent variable, y is the dependent variable, and a is
    % the model parameter.
    
    a_min = a_range(1);
    a_max = a_range(2);
    max_iters = 10;
    for j = 1:max_iters
        a = linspace(a_min,a_max,100)'; % search space for parameter
        loss = sum((y-model(a,x)).^2,2);
        [~,opt_id] = min(loss);
        a_est = a(opt_id);
        a_min = a(max(opt_id-1,1));
        a_max = a(min(opt_id+1,numel(a)));
    end
end
