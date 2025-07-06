function xyb = GenerateRandomConstellation(n, min_sep_frac, contrast, scale)
    % Randomly generates the coordinates and brightnesses for a 
    % constellation consisting of n point sources located within a disk of 
    % radius 0.5 wherein each source has at least one nearest neighbor that 
    % is min_sep away. The contrast sets the ratio of the brightest source
    % to the dimmest one. The center of intensity (center-of-mass) of the 
    % constellation lies at the origin (0,0). The brightnesses are given in
    % relative scale such that they sum to one.
    % ------------------
    % INPUTS:
    % ------------------
    % n              : number of point sources in the constellation
    % min_sep_frac   : minimum separation distance in constellation. Must lie in the interval (0,1)
    % contrast       : ratio of minimum brightness to the maximum
    %                  brightness (contrast = min_b/max_b) which lies in the interval (0,1]
    % scale          : scales the coordinates of the sources
    % ------------------
    % OUTPUTS
    % ------------------
    % xyb             : nx3 matrix containing the xy and brightness coordinates of the
    %                  constellation. 
    %                       xyb(:,1) is the x coordinate
    %                       xyb(:,2) is the y coordinate 
    %                       xyb(:,3) is the b coordinates (sums to 1)
                            
    
    % half-unit support disk
    D = 1;          % Diameter of support disk
    R = D/2;        % Radius of support disk
    
    % assert that more than one source was requested
    assert(n>=1);

    % assert that the contrast is valid
    assert(0<contrast && contrast<=1)

    %packing fractions for circle-within-a-circle up to 20 circles
    % https://en.wikipedia.org/wiki/Circle_packing_in_a_circle
    %circ_pack_frac = [1,0.5000,0.6466,0.6864,0.6854,0.6666,0.7777,0.7328,0.6895,0.6878,0.7148,0.7392,0.7245,0.7474,0.7339,0.7512,0.7403,0.7609,0.8034,0.7623];

    % no possible samples if the area of the area from the optimal packing fraction
    % is less than the area that the sources may be.
    %assert(n*(min_sep_frac/2).^2 <= circ_pack_frac(n) *  (R+min_sep_frac/2)^2); 
    
    % uniformly sample a point anywhere within the disk  of radius R-min_sep 
    % (ensures that the second point is guaranteed to lie within the disk
    % of radius R)
    p1 = sampleDiskUniform(1,R-min_sep_frac); 
    
    % add p1 to the list of point coordinates
    xy(1,:) = p1;    
        
    % variablility in deviation from min_sep
    epsilon = min_sep_frac/100;

    % generate remaining samples
    for k = 2:n
        % check if all the points are within the min separation criteria,
        % otherwise regenerate the kth point
        while size(xy,1) < k || ~all(pdist(xy(1:k,:)) >= min_sep_frac - (epsilon/2))
            
            % sample a point on a fuzzy ring of width epsilon with radius min_sep
            [rkx,rky] = pol2cart(2*pi*rand(1), epsilon*(rand(1)-.5) + min_sep_frac); 
            rk = [rkx,rky];
            
            % randomly pick a point to add the sampled ring point coordinates
            j = randi(k-1);
            pj = xy(j,:);   
            
            % generate the new point pk
            pk = pj + rk;
            
            % add pk to the list of point coordinates after the random
            % index j
            xy(k,:) = pk;  

        end
    end

    % sample the brightnesses
    b = unifrnd(contrast,1,[n-2,1]);
    b = [contrast; b; 1];
    b = b(randperm(n));
    b = b/sum(b);
    
    % realign the centroid
    xy = xy - sum(b.*xy,1);
    
    % check if the scene still falls inside the FOV. Otherwise rerun
    % the function
    %{
    if any( sum(xy.^2,2) > R^2)
        xyb = GenerateRandomConstellation(n, min_sep_frac, contrast, scale);
    else
        xyb = [scale*xy,b];
    end
    %}
    xyb = [scale*xy,b];

    % sort coordinates in ascending order of brightness
    [~,idx] = sort(xyb(:,3),1,"ascend");
    xyb = xyb(idx,:);   
end

function xy = sampleDiskUniform(n,R)
   % generates n samples uniformly over the unit disk
   r = R*sqrt(rand(n,1));
   th = 2*pi*rand(n,1);
   
   [x,y] = pol2cart(th,r);
   xy = [x,y];
end