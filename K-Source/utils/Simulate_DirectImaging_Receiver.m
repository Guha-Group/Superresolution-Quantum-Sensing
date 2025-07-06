function [xyb_est,out_data] = Simulate_DirectImaging_Receiver(xyb,M,N,sigma,visualization_flag)
     
    % extract parameters
    num_sources = size(xyb,1);
    xy = xyb(:,1:2);
    b = xyb(:,3);

    %% CALIBRATION STAGE
    % set the calibration scene (equal brightness)
    b0 = ones(num_sources,1)/num_sources;
    xyb0 = [xy,b0];

    % direct imaging measurement to get centroid
    xy_samples = DirectImagingMeasurement(xyb0,sigma,M);

    % estimate the source positions
    xy_est = DILocalizeSources(xy,xy_samples,sigma);
    
    %% SENSING STAGE

    % direct imaging measurement to get brightness pre-estimates
    xy_samples = DirectImagingMeasurement(xyb,sigma,N);

    % estimate brightness from direct imaging measurement
    b_est = EstimateBrightnessDI(xy_samples,xy_est,sigma);

    % collect all estimates
    xyb_est = [xy_est,b_est];

    % return the collection of estimates (sorted by brightness from dimmest to brightest)
    xyb_est = [xy_est, b_est];
    %[~,id] = sort(b_est,'ascend');
    %xyb_est = xyb_est(id,:);

    out_data = [];

    %--------------------------------------------------%
    if visualization_flag
        figure
        T = tiledlayout(1,3,"TileSpacing","compact","Padding","compact");
        title(T,'Direct Imaging Receiver','interpreter','latex')

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