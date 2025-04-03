function [x0_est,mom2_est,M1_mom2] = DirectImagingParameterPreEstimates(x,sigma)
     % This function attempts to robustly determine estimate of the
     % centroid and the second moment of the distribution (related to the
     % source pair separation distance) given a list of direct imaging
     % measurement samples x. We wish to enforce the requirement that
     % mom2_est>1 which corresponds to a situation where there is
     % definitively a non-zero separation between two point sources.

     % Attempt 1: 
     % Try the basic moment approximations (guaranteed to work for large enough numbers of samples)
     x0_est = mean(x);
     mom2_est = mean(((x-x0_est)/sigma).^2);
     M1_mom2 = numel(x);

     % check if the whitened second moment is less than 1 
     if mom2_est<1
         
        % Attempt 2: 
        % Try solving system of nonlinear equations for solving simultaneous
        % maximim-likelihod estimator of midpoint and separation
        DI_MLE = @(mle) [sum(x + mle(2)* tanh( (x-mle(1)) * mle(2) / sigma^2) ) / M1_mom2 - mle(1);...
                         sum( tanh( (x-mle(1)) * mle(2) / sigma^2) .*(x-mle(1))) - mle(2)];
        options = optimset('Display','off');
        mle_out = fsolve(DI_MLE,[0,sigma*.5],options);
        x0_mle = mle_out(1);    
        s_mle = mle_out(2);

        % check if the MLE solution is valid 
        if s_mle<1e-7 % smallest value possible such that floating point precision doesn't snap to zero when squaring
            
            
            % Attempt 3:
            % try peeling off samples near the midpoint until the second moment is
            % larger that 1 (biases the samples to concentrate away from
            % the first moment)
            
            % estimate midpoint
            x0_est = mean(x); 
            
            % sort the samples in ascending order based on their absolute 
            % distance to the midpoint
            [~,dist_to_midpoint] = sort(abs(x-x0_est)); 
            x_reduced = x(dist_to_midpoint);

            % while the second moment is less than 1 and we have samples
            % left...
            while (mom2_est) < 1 && (numel(x_reduced) > 0)
                % peel off a sample
                x_reduced = x_reduced(2:end);
                % compute the second moment on the reduced sample set
                mom2_est = mean(((x_reduced-x0_est)/sigma).^2);
            end
            
            % check if the whitened second moment is less than one
            if mom2_est < 1

                % Attempt 4: 
                % At this point we have no other options but to set the
                % second moment to be some small perturbation greater than
                % one.
                mom2_est = 1+ 1e-6;
            else
                M1_mom2 = numel(x_reduced);
            end
        else
            % assign the outputs
            x0_est = x0_mle;
            mom2_est = 1 + (s_mle/sigma)^2;
        end
     end
end