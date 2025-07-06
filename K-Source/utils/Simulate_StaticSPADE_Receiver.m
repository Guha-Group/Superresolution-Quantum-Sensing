function [xyb_est,out_data] = Simulate_StaticSPADE_Receiver(xyb,M,N,sigma,splitting_ratio,nm,visualize_flag)
    % Simulates a Multi-stage SPADE-enhanced receiver for performing
    % sub-diffraction color center sensing. 
    % DESCRIPTION: The receiver executes a
    % calibration stage and a sensing stage in sequence. Each of these two
    % stages are further decomposed into a direct imaging measurement and
    % an a SPADE measurement. The calibration stage acts on the
    % constellation configuration driven to equal brightness. Here is where
    % the receiver estimates the locations of the emitters. In the sensing
    % stage the receiver estimates the brightnesses of the sources
    % conditioned on the estimated locations.
    %
    %%%%%%%%%%%%%%%%%%%
    %%%%% INPUTS %%%%%%
    %%%%%%%%%%%%%%%%%%%
    % xyb               : [K,3] scene array with each row contianing the xy
    %                     coordinate and relative brightness of each emitter
    % M                 : number of calibration stage photons
    % N                 : number of sensing stage photons
    % sigma             : Rayleigh unit
    % splitting_ratio   : fraction of photons to direct imaging measurements
    % nm                : [M,2] HG mode indices
    % visualize_flag    : boolean value indicating whether results are to
    %                     be displayed or not.
    %%%%%%%%%%%%%%%%%%%
    %%%% OUTPUTS %%%%%%
    %%%%%%%%%%%%%%%%%%%
    % xyb_est           : [K,3] estimated scene parameter array
    % out_data          : a structure containing output metadata

    
    % extract parameters
    num_sources = size(xyb,1);
    xy = xyb(:,1:2);
    b = xyb(:,3);

    %% CALIBRATION STAGE
    % set the calibration scene (equal brightness)
    b0 = ones(num_sources,1)/num_sources;
    xyb0 = [xy,b0];
    
    % allocate photons for calibration stage
    M1 = ceil(splitting_ratio*M);
    M2 = M-M1;

    % direct imaging measurement to get centroid
    xy_samples = DirectImagingMeasurement(xyb0,sigma,M1);

    % centroid estimate
    centroid_est = mean(xy_samples,1);

    % change pointing 
    xy = xy - centroid_est;
    xyb0 = [xy,b0];

    % SPADE measurement
    rot_angles = [0,pi/4];
    nm_samples = HGSPADEMeasurement(xyb0,nm,M2,sigma,rot_angles);

    % estimate the source positions from the spade measurement
    xy_est = HGLocalizeSources(xy,nm,nm_samples,sigma,rot_angles,visualize_flag);    
   
    %% SENSING STAGE

    % set the sensing scene (unequal brightness)
    xyb1 = [xy,b];
    
    % allocate photons for sensing stage
    N1 = ceil(splitting_ratio*N);
    N2 = N-N1;

    % direct imaging measurement to get brightness pre-estimates
    xy_samples = DirectImagingMeasurement(xyb1,sigma,N1);

    % estimate brightness from direct imaging measurement
    b_pre_est = EstimateBrightnessDI(xy_samples,xy_est,sigma);

    % formulate preliminary scene estimate
    xyb_est = [xy_est,b_pre_est];
    
    % formulate YKL measurement
    [YKL_samples,YKL,Psi_est,YKL_GaussRepn] = YKLMeasurement(xyb1,xyb_est,sigma,N2,visualize_flag);

    % estimate the source brightnesses 
    b_est = EstimateBrightnessYKL(YKL_samples,YKL,Psi_est);
    
    % revert the estimates to their original coordinate system
    xy_est = xy_est + centroid_est;
    
    % return the collection of estimates (sorted by brightness from dimmest to brightest)
    xyb_est = [xy_est, b_est];
    %[~,id] = sort(b_est,'ascend');
    %xyb_est = xyb_est(id,:);
    
    % store any desired meta-data
    out_data.YKL_GaussRepn = YKL_GaussRepn;


    %--------------------------------------------------%
    % visualize estimates
    if visualize_flag
        figure
        T = tiledlayout(1,3,"TileSpacing","compact","Padding","compact");
        title(T,'SPADE Receiver','interpreter','latex')
        
        % show localization estimates
        nexttile(1)
        hold on
        scatter(xyb(:,1),xyb(:,2),'k','filled')
        scatter(xyb_est(:,1),xyb_est(:,2),'r','d')
        hold off
        axis square; box on; grid on;
        xticks(-1:.5:1); yticks(-1:.5:1);
        xlim([-1,1]); ylim([-1,1]);
        xlabel('$x/\sigma$','interpreter','latex')
        ylabel('$y/\sigma$','interpreter','latex')
        title('Location Estimates','interpreter','latex')
        legend({'Ground Truth','Estimate'},'interpreter','latex')
    
        % show brightness estimates
        nexttile(2)
        bar_chart = bar(1:num_sources,[xyb(:,3),xyb_est(:,3)],'FaceAlpha',.5);
        bar_chart(1).FaceColor = [0,0,0];
        bar_chart(2).FaceColor = [1,0,0];
        axis square; box on;
        xlabel('Source Index','interpreter','latex')
        ylabel('Relative Brightness','interpreter','latex')
        title('Brightness Estimates','interpreter','latex')
        legend({'Ground Truth','Estimate'},'interpreter','latex')
        
        % show reconstructed scene
        nexttile(3)
        hold on
        scatter(xyb(:,1),xyb(:,2),30*size(xyb,1)*xyb(:,3),'k','filled')
        scatter(xyb_est(:,1),xyb_est(:,2),30*size(xyb_est,1)*xyb_est(:,3),'r','d')
        hold off
        axis square; box on; grid on;
        xticks(-1:.5:1); yticks(-1:.5:1);
        xlim([-1,1]); ylim([-1,1]);
        xlabel('$x/\sigma$','interpreter','latex')
        ylabel('$y/\sigma$','interpreter','latex')
        title('Scene Reconstruction','interpreter','latex')
        legend({'Ground Truth','Estimate'},'interpreter','latex')
    end
    %--------------------------------------------------%


end