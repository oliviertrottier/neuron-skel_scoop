function [Tree, Soma_PixelInd, image_bin] = skel_scoop(image_input, varargin)
% Function to skeletonize a neuron by scooping voxels. The algorithm is
% based the voxel scooping algorithm developped by
% Rodriguez, Alfredo, Douglas B. Ehlenberger, Patrick R. Hof, and Susan L. Wearne. “Three-Dimensional Neuron Tracing by Voxel Scooping.” Journal of Neuroscience Methods 184, no. 1 (October 30, 2009): 169–75. https://doi.org/10.1016/j.jneumeth.2009.07.021.

% Input
% image_input = 2D array representing the image.

% Output
% Tree = structure containing all branches of the neuronal tree. Each entry
%       in the structure corresponds to one branch.
%       The fields of Tree are:
%       1) PointsPos = Positions of the points forming the branch.
%       2) PointsDiameter = Approximate diameter of the branch at each point position.
%       3) Length = Number of segments that form the branch.
%       4) Length_dim = Total dimensionfull length of the branch (given in units of pixels).
%       5) ParentID = ID of the parent branch.
%       6) SiblingID = ID of the sibling (or sister) branch.
%       7) ChildrenID = IDs of the children (or daughter) branches.
%       8) Depth = Depth of the branch in the tree structure. The root branches have a depth of 0.
% Soma_PixelInd = 2D index of the Soma pixel.
% image_bin = binarized version of the input image used for scooping.
%% Parse optional parameters
p = inputParser;
addParameter(p, 'SomaPos', [0, 0]); % Fix the position of the Soma in the output structure.
addParameter(p, 'Seed', 'Auto'); % Index position of the seed.
addParameter(p, 'Close', false); % Perform morphological closure before tracing. Useful to remove holes.
addParameter(p, 'Plots', nargout==0); % Plot intermediate plots detailing the skeletonization steps.
addParameter(p, 'PlotIterations', false); % Plot the cluster positions found in each iteration of the voxel scooping step.
addParameter(p, 'Waitbar', true); % Display waitbar.
addParameter(p, 'L_min', 'local_diameter'); % Min branch length (in pixels). By default, it uses the local diameter of the branch as the threshold.
addParameter(p, 'Clustersize_min', 1); % Minimum number of pixels that a cluster must have.
addParameter(p, 'Pixel_mass', 1); % Pixel mass used to measure the center-of-mass of each cluster.
addParameter(p, 'Pixel_size', 1); % Pixel size used to rescale the skeleton nodes position.
addParameter(p, 'Branch_nodes_distance', 1); % Distance between two nodes of a branch in pixels (dimensionless).
addParameter(p, 'Binarization_method', 'local_threshold'); % Threshold for local thresholding (expressed as a proportion of the local mean).
addParameter(p, 'Binarize_only', false); % Return only the binarized image.
addParameter(p, 'Mean_filter_thresh', 1.225); % Threshold for local mean thresholding (expressed as a proportion of the local mean).
addParameter(p, 'Soma_mask', 0); % Calculate a mask that encloses the soma. (Not implemented).
parse(p, varargin{:});
options = p.Results;
%% Find the location of the seed (Soma).
image_size = size(image_input);
image_intmax = double(intmax(class(image_input)));
if isa(options.Seed, 'char') && strcmp(options.Seed, 'Auto')
    if islogical(image_input)
        % If the image is a logical, use the gray-weigthed distance
        % transform to find the pixel that is the furthest from background pixels.
        image_GWDT = bwdist(~image_input);
        [~, seed_ind] = max(image_GWDT(:));
        seed_ind = [mod(seed_ind-1, image_size(1))+1, ceil(seed_ind/image_size(1))];
        
    else
        % If the image is not binary, binarize the image to retain bright
        % pixels. Mean-filter the binary image and choose the pixel with the highest
        % intensity.
        
        % Binarize the image globally to retain a minimum of 5% bright pixels.
        saturation_levels = stretchlim(image_input,0.05);
        bright_thres = min(0.99,saturation_levels(2));
        image_bin_bright = imbinarize(image_input,bright_thres);
        
        % Fill holes in the bright binary image. This is sometimes
        % necessary when soma centers have lower pixel intensities in their
        % center.
        %image_bin_bright = imfill(image_bin_bright,'holes');
        
        % Calculate the gray-weigthed image transform.
        image_GWDT = bwdist(~image_bin_bright);
        
        % Theshold the GWDT image to find a mask that identifies potential soma positions. The
        % center of soma are usually deeper than 10 pixels.
        Soma_depth_min = 10;
        soma_mask = image_GWDT > Soma_depth_min;
        
        % Find all the connected components in the soma_mask.
        Soma_CC = bwconncomp(soma_mask);
        
        % Threshold the soma candidate region using their Area.
        CC_N_pixels = cellfun(@numel,Soma_CC.PixelIdxList);
        Soma_CC.PixelIdxList = Soma_CC.PixelIdxList(CC_N_pixels > 500);
        Soma_CC.NumObjects = numel(Soma_CC.PixelIdxList);
        N_soma = Soma_CC.NumObjects;
        
        if N_soma < 0
            error('A soma could not be found');
        elseif N_soma > 1
            % If there are
            % more than one connected component, choose the connected component
            % whose centroid is closest to the center of the image.
            disp('More than 1 seed was found. The most centered one is chosen.')
            
            CC_props = regionprops(Soma_CC);
            Image_center = round(image_size/2);
            for i=1:N_soma
                CC_props(i).DistanceToSeed = sqrt(sum((CC_props(i).Centroid - Image_center).^2,2));
            end
            [~,Soma_ind] = min([CC_props.DistanceToSeed]);
            Soma_Region_pixel_ind = Soma_CC.PixelIdxList{Soma_ind};
            [~, seed_ind] = max(image_GWDT(Soma_Region_pixel_ind));
            seed_ind = Soma_Region_pixel_ind(seed_ind);
        else
            [~, seed_ind] = max(image_GWDT(:));
        end
        
        % Change the seed linear index to a 2D index.
        seed_ind = [mod(seed_ind-1, image_size(1))+1, ceil(seed_ind/image_size(1))];
    end
    
    % Plot the image with the seed mask.
    if options.Plots
        seed_mask = false(img_size);
        seed_mask(seed_ind(1),seed_ind(2)) = true;
        figure;imshowpair(image_input,seed_mask);
    end
    
else
    seed_ind = options.Seed;
end
Soma_PixelInd = seed_ind;
%% Binarize the input image (if not already binarized).
if ~islogical(image_input)
    if options.Waitbar
        wb = waitbar(0,'Binarizing the image.');
    end
    
    % Create an image where zero-intensity pixels are set to background intensity.
    % This is done to avoid the effects of cropped regions
    % (intensity=0) on the calculation of the local mean.
    % The background intensity is found by finding the most frequent intensity
    % value that is below the otsu threshold of the image.
    nonzero_pixel_intensities = image_input(image_input>0);
    bin_edges = 0:double(max(nonzero_pixel_intensities));
    pixel_counts = histc(double(nonzero_pixel_intensities),bin_edges);
    otsu_thresh = round(otsuthresh(pixel_counts)*numel(pixel_counts));
    [~,max_ind] = max(pixel_counts(1:otsu_thresh));
    %mean_pixel_intensity = round(mean(nonzero_pixel_intensities));
    image_input_mod = image_input;
    image_input_mod(image_input_mod == 0) = max_ind;
    
    % Adjust the contrast.
    saturation_levels = stretchlim(image_input_mod,[0.01 0.975]);
    image_adj = imadjust(image_input_mod,saturation_levels,saturation_levels);
    
    % Filter the image with a Wiener filter to remove noise.
    im_filtered = wiener2(image_adj,3*ones(1,2));
    
    switch options.Binarization_method
        case 'local_threshold'
            %% Computer a local threshold using the mean as the threshold statistic.
            filter_size = 51;
            %filter_size = 31;
            mean_filter = ones(filter_size)/filter_size^2;
            %mean_thresh = 1.225;
            %mean_thresh = 1.3;
            mean_thresh = options.Mean_filter_thresh;
            
            padSize = (filter_size-1)/2;
            im_padded = padarray(im_filtered,padSize*ones(1,2),'replicate','both');
            %im_padded = padarray(image_adj,padSize*ones(1,2),'replicate','both');
            
            threshold_matrix = mean_thresh/double(intmax(class(image_input)))*conv2(double(im_padded),mean_filter,'valid');
            
            %local_threshold2 = adaptthresh(im_filtered,1.6 - mean_thresh,'NeighborhoodSize',filter_size,'Statistic','mean','ForegroundPolarity','bright');
            %local_threshold = adaptthresh(im_filtered,0.5,'NeighborhoodSize',filter_size,'Statistic','mean');
            %local_threshold = adaptthresh(im_filtered,0.1,'NeighborhoodSize',filter_size,'Statistic','median');
            
            % Plot the image binarized with the local threshold.
            %image_bin = imbinarize(image_input,local_threshold);
            %im_adj = imadjust(image_input);
            %figure;imshowpair(im_adj,image_bin);
        case 'global_threshold_noise_fit'
            %% Fit a gaussian to the noise and use the fit to define a global threshold.
            %pixel_intensities = im_filtered(:);
            pixel_intensities = image_input(image_input>0);
            bin_edges = (0:double(max(pixel_intensities)+1))';
            bin_centers = bin_edges(1:end-1) + (bin_edges(2)-bin_edges(1))/2;
            pixel_counts = histc(double(pixel_intensities),bin_edges);
            pixel_counts = pixel_counts(1:end-1);
            [max_counts,peak_ind] = max(pixel_counts);
            peak_pos = bin_centers(peak_ind);
            
            gaussian_fit_model = @(beta,x) beta(1).*exp(-(x - peak_pos).^2./(2*beta(2)^2))./sqrt(2*pi)/beta(2);
            fit_weights = [ones(peak_ind,1); 0*ones(numel(bin_centers)-peak_ind,1)];
            %beta_fit = nlinfit(bin_centers(1:peak_ind),pixel_counts(1:peak_ind),gaussian_fit_model,[max_counts,10],'Weights',fit_weights);
            
            % Fit the left part of the peak.
            gaussian_fit_beta = nlinfit(bin_centers(1:peak_ind),pixel_counts(1:peak_ind),gaussian_fit_model,[max_counts,10]);
            
            % Fit the right part of the peak.
            gaussian_fit_beta = nlinfit(bin_centers(peak_ind:end),pixel_counts(peak_ind:end),gaussian_fit_model,[max_counts,10]);
            
            gaussian_noise_threshold = round(peak_pos + 3*gaussian_fit_beta(2));
            im_bin = imbinarize(image_input,gaussian_noise_threshold/image_intmax);
            
            % Fit an exponential to the pixel counts above threshold.
            exp_fit_model = @(beta,x) beta(1).*exp(-x/beta(2));
            bin_above_threshold = bin_centers > gaussian_noise_threshold;
            exp_fit_bin_centers = bin_centers(bin_above_threshold);
            exp_fit_beta = nlinfit(exp_fit_bin_centers,pixel_counts(bin_above_threshold),exp_fit_model,[max_counts,10]);
            
            % Define the threshold matrix used to threshold the image.
            threshold_matrix = gaussian_noise_threshold/double(intmax(class(image_input)))*ones(size(im_bin));
            
            if options.Plots
                % Plot the histogram.
                figure;hold on;
                plot(bin_centers,pixel_counts,'x','DisplayName','Histogram');
                bin_centers_fit = 0:.1:bin_centers(end);
                
                % Plot gaussian fit around peak.
                plot(bin_centers_fit,gaussian_fit_model(gaussian_fit_beta,bin_centers_fit),'-','DisplayName','Gaussian Fit')
                a = gca;
                
                % Plot threshold.
                plot(gaussian_noise_threshold*ones(1,2),a.YLim,'-','DisplayName','Threshold');
                
                % Plot the exponential fit for the pixel counts above threshold.
                plot(exp_fit_bin_centers,exp_fit_model(exp_fit_beta,exp_fit_bin_centers),'--','DisplayName','Above Threshold Exp. Fit');
                
                legend;
                
                % Plot an overlay of the binary image with the contrast
                % adjusted image.
                figure;imshowpair(image_adj,im_bin);a=gca;tightax;a.Visible='on';
            end
        case 'contour_contraction'
            %% Find intensity contours in the image and use them to define the neuron mask.
            conn_comps = bwconncomp(im_filtered > 120);
            im_biggest_component = (min(im_filtered(:)))*ones(size(im_filtered),class(im_filtered));
            im_biggest_component(conn_comps.PixelIdxList{1}) = im_filtered(conn_comps.PixelIdxList{1});
            
            % Find the threshold where the main connected component starts
            % breaking apart.
            [N_conn_comps,levels] = calculate_conn_comps(im_biggest_component,50);
            threshold_ind = find(diff(N_conn_comps) >= 4,1);
            %threshold_ind = max(1,find(N_conn_comps > 1,1)-1);
            global_threshold = levels(threshold_ind);
            
            im_bin2 = im_biggest_component > levels(threshold_ind+1);
            im_conn_comps = bwlabel(im_bin2);
            %figure;imagesc(im_conn_comps);
            %figure;imshowpair(im_biggest_component,im_bin2);
            
            threshold_matrix = global_threshold/double(intmax(class(image_input)))*ones(size(im_filtered));
        case 'edge_detection'
            %% Binarize the image by first finding edges and fill the interior region circumscribed by the edges.
            %edge_threshold = 1e-4;
            %im_edges = edge(im_filtered,edge_threshold);
            %figure;imshowpair(im_filtered,im_edges);
            [im_edges_filled,im_edges] = im_gradient_fill(im_filtered);
            
            threshold_matrix = ones(size(im_filtered));
            threshold_matrix(im_edges_filled) = 0;
            
            %figure;imshowpair(im_filtered,im_edges_filled);
            %im_edges_filled_closed = bwmorph(im_edges_filled,'close');
            %figure;imshowpair(im_filtered,im_edges_filled_closed);
        case 'global_otsu'
            %% Define the neuron mask as all pixels whose intensity is above a global threshold defined by otsu's method.
            max_pixel_val = double(max(im_filtered(:)));
            bin_edges = 0:1:max_pixel_val;
            bincounts = histc(im_filtered(:),bin_edges);
            otsu_thresh = otsuthresh(bincounts);
            global_threshold = otsu_thresh * max_pixel_val;
            threshold_matrix = global_threshold/double(intmax(class(image_input)))*ones(size(im_filtered));
            
        case 'global_multi_thresh'
            %% Define the neuron mask as all pixels whose intensity is above a global threshold defined by otsu's multi-threshold method.
            multi_thresholds = multithresh(im_filtered,2);
            im_quantized = imquantize(im_filtered,multi_thresholds);
            threshold_matrix = double(~(im_quantized==3));
    end
    
    % Binarize the image with a high threshold. This is useful to keep clusters of high
    % intensity pixels that may be removed by local thresholding.
    saturation_levels = stretchlim(image_input_mod,[0.01 0.95]);
    bright_threshold = min(0.99,saturation_levels(2));
    image_bin_bright = imbinarize(image_input_mod,bright_threshold);
    
    % Compute an approximate size of the soma using the seed position.
    Soma_size_avg = find_branch_thickness(seed_ind,image_bin_bright,'Statistic','mean');
    Soma_radius_avg = ceil(Soma_size_avg/2);
    
    % Recompute the threshold in the region close to the seed. This region
    % is usually bright and local thresholding can remove clusters of
    % bright pixels, if they are big enough.
    Soma_region_radius = Soma_radius_avg;
    Soma_row_ind = seed_ind(1) + (-Soma_region_radius:Soma_region_radius);
    Soma_col_ind = seed_ind(2) + (-Soma_region_radius:Soma_region_radius);
    threshold_matrix(Soma_row_ind,Soma_col_ind) = min(threshold_matrix(Soma_row_ind,Soma_col_ind),bright_threshold);
    
    % Binarize image with the threshold matrix.
    image_bin = imbinarize(image_input_mod,threshold_matrix);
    
    % Fill holes in the soma region.
    image_bin(Soma_row_ind,Soma_col_ind) = imfill(image_bin(Soma_row_ind,Soma_col_ind),'holes');
    
    % Remove isolated pixels.
    image_bin = bwmorph(image_bin,'clean');
    
    % Find all connected components in the binary image.
    CC = bwconncomp(image_bin);
    Components_pixels_ind = CC.PixelIdxList';
    
    % Remove pixels that are part of small clusters.
    Components_N_pixels = cellfun(@numel,Components_pixels_ind);
    Small_clusters = Components_pixels_ind(Components_N_pixels < max(Components_N_pixels)/2);
    Small_clusters = cell2mat(Small_clusters);
    image_bin(Small_clusters) = 0;
    
    if options.Waitbar
        waitbar(1,wb);
        delete(wb);
    end
else
    image_bin = image_input;
end

% Throw error if the seed is not found in the binary image.
if ~image_bin(seed_ind(1),seed_ind(2))
    error('The binary image does not contain the seed');
end

% Close the binary image to remove holes.
if options.Close
    image_bin = bwmorph(image_bin, 'close');
end

% Plot the final binary image overlaid on top of the contrast-adjusted image.
if options.Plots
    %figure;imshow(image_bin)
    figure;imshowpair(image_adj,image_bin);a=gca;tightax;a.Visible='on';
end

% End early if only a binarized image is necessary.
if options.Binarize_only
    Tree = [];
    Soma_PixelInd = [];
    return;
end
%% Find the soma mask.
if options.Soma_mask
    error('(under development)');
    [soma_candidate_row,soma_candidate_col] = ndgrid(-50:50,-50:50);
    soma_candidate_pixel_ind = [soma_candidate_row(:)+seed_ind(1), soma_candidate_col(:)+seed_ind(2)];
    soma_candidate_pixel_lin_ind = soma_candidate_pixel_ind(:,1) + (soma_candidate_pixel_ind(:,2)-1)*image_size(1);
    
    soma_candidate_pixel_ind = soma_candidate_pixel_ind(image_bin(soma_candidate_pixel_lin_ind),:);
    soma_candidate_pixel_lin_ind = soma_candidate_pixel_lin_ind(image_bin(soma_candidate_pixel_lin_ind));
    soma_candidate_thickness = find_branch_thickness(soma_candidate_pixel_ind, image_bin);
    
    image_bin_soma = image_bin((-50:50) + seed_ind(1),(-50:50) + seed_ind(2));
    image_thickness = zeros(image_size);
    image_thickness(soma_candidate_pixel_lin_ind) = soma_candidate_thickness;
    soma_mask = image_thickness >= max(image_thickness(:))/4;
    
    image_edges = edge(im_filtered);
    figure;imagesc(image_thickness);
    figure;imshowpair(soma_mask,image_bin);
end
%% Define the pixel mass.
% The pixel mass is used to find the clusters' center of mass. The center
% of mass orients the next step of the scooping step.
use_pixel_mass = options.Pixel_mass ~= 1;

if isnumeric(options.Pixel_mass) && numel(options.Pixel_mass) > 1
    options.Pixel_mass = double(options.Pixel_mass);
elseif strcmp(options.Pixel_mass,'GWDT')
    image_GWDT = bwdist(~image_bin);
    options.Pixel_mass = image_GWDT;
elseif strcmp(options.Pixel_mass,'GWDT_log')
    image_GWDT = bwdist(~image_bin);
    options.Pixel_mass = log(image_GWDT);
elseif strcmp(options.Pixel_mass,'Intensity')
    options.Pixel_mass = double(image_input);
end
%% Initialization
N_clusters_per_iter = 100;
N_iterations_init = 1000; % The number of iterations is free to increase and is mainly used to set initialize sizes of arrays.
Cluster_size_min = options.Clustersize_min; % Minimum number of pixels that a cluster must have.
neigborhood_radius = 1; % Size of a cluster's neighborhood.
search_radius = 5; % Search radius (in pixel) used to find the next cluster pixels.

N_dims = ndims(image_bin);
im_size = size(image_bin);
N_pixels = nnz(image_bin);

% Keep track of pixels that remain to be visited (scooped).
is_pixel_notvisited = image_bin;
pixels_clusterID = zeros(size(image_bin));

Clusters_pixel_inds = cell(N_iterations_init, N_clusters_per_iter);
Clusters_size = nan(N_iterations_init, N_clusters_per_iter);
Clusters_position = nan(N_iterations_init, N_clusters_per_iter, N_dims);
Clusters_ID = nan(N_iterations_init, N_clusters_per_iter);
Clusters_parentID = nan(N_iterations_init, N_clusters_per_iter);
N_clusters = zeros(N_iterations_init, 1);

% Initialize arrays with seed.
is_pixel_notvisited(seed_ind(1), seed_ind(2)) = false;
pixels_clusterID(seed_ind(1), seed_ind(2)) = 1;
Clusters_pixel_inds{1, 1} = seed_ind;
Clusters_size(1, 1) = sqrt(2);
Clusters_position(1, 1, :) = seed_ind - 0.5*ones(1, N_dims); % Initialize at the center of the seed pixel.
Clusters_ID(1, 1) = 1;
Clusters_parentID(1, 1) = 0;
N_clusters(1) = 1;
Cluster_ID = 1;

% Start at the second iteration (iteration = 1 is the initialization).
i = 2;
%% Skeletonization
if options.Waitbar
    wb = waitbar(0, 'Skeletonizing neuron...');
end

while N_clusters(i-1) > 0
    if options.Waitbar
        waitbar(1-nnz(is_pixel_notvisited)/N_pixels, wb);
    end
    
    if options.PlotIterations
        pixels_clusterID_previous = pixels_clusterID;
    end
    
    % Loop over clusters found in the previous iteration to find the new
    % clusters.
    N_new_clusters = 0;
    for j = 1:N_clusters(i-1)
        % Get the pixel indices corresponding to cluster C(i-1,j)
        cluster_pixels_ind = Clusters_pixel_inds{i-1, j};
        
        % Define the range to search for unvisited pixels in the
        % neighborhood of the cluster.
        rowcol_minmax = [min(cluster_pixels_ind, [], 1); max(cluster_pixels_ind, [], 1)];
        rowcol_minmax = rowcol_minmax + search_radius*[-1, -1;1, 1];
        rowcol_minmax(1, :) = max(rowcol_minmax(1, :), [1, 1]);
        rowcol_minmax(2, :) = min(rowcol_minmax(2, :), im_size);
        row_inds = rowcol_minmax(1, 1):rowcol_minmax(2, 1);
        col_inds = rowcol_minmax(1, 2):rowcol_minmax(2, 2);
        
        % Find unvisited pixels in a small region near the cluster.
        is_region_pixel_notvisited = is_pixel_notvisited(row_inds, col_inds);
        [region_row_ind, region_col_ind] = find(is_region_pixel_notvisited);
        Region_pixels_inds = [region_row_ind, region_col_ind];
        Region_pixels_global_inds = bsxfun(@plus, Region_pixels_inds, rowcol_minmax(1, :)-1);
        
        % Find unvisited pixels that are neighbours of the cluster.
        [Neighbour_inds, is_region_pixel_in_neighborhood] = find_cluster_neigborhood(cluster_pixels_ind, Region_pixels_global_inds, neigborhood_radius);
        
        %         image_cluster = false(im_size);
        %         Neighbour_lin_inds = cluster_pixels_ind(:,1) + (cluster_pixels_ind(:,2) - 1)*im_size(1);
        %         image_cluster(Neighbour_lin_inds) = true;
        %         figure;imshowpair(image_cluster,image_bin);
        %
        %         image_neighbourhood = false(im_size);
        %         Neighbour_lin_inds = Neighbour_inds(:,1) + (Neighbour_inds(:,2) - 1)*im_size(1);
        %         image_neighbourhood(Neighbour_lin_inds) = true;
        %         figure;imshowpair(image_neighbourhood,image_bin);
        %
        %         image_region = false(im_size);
        %         Region_pixels_lin_inds = Region_pixels_global_inds(:,1) + (Region_pixels_global_inds(:,2) - 1)*im_size(1);
        %         image_region(Region_pixels_lin_inds) = true;
        %         figure;imshowpair(image_region,image_bin);
        
        % Move to the next cluster if no neighbours are found.
        if isempty(Neighbour_inds)
            continue;
        end
        
        % Find the connected clusters among the cluster's neighbours.
        [Neighbours_cluster_ID, N_connected_clusters, Connected_clusters_size] = find_connected_clusters(Neighbour_inds);
        
        % Remove neighbour pixels that belong to clusters whose size is below
        % the minimal size.
        Connected_cluster_IDs = 1:N_connected_clusters;
        is_connected_cluster_small = Connected_clusters_size < Cluster_size_min;
        Small_clusters_ID = Connected_cluster_IDs(is_connected_cluster_small);
        if ~isempty(Small_clusters_ID)
            is_Neighbours_in_small_cluster = builtin('_ismemberhelper', Neighbours_cluster_ID, Small_clusters_ID);
            Neighbours_cluster_ID(is_Neighbours_in_small_cluster) = [];
            Neighbour_inds(is_Neighbours_in_small_cluster, :) = [];
            
            Connected_cluster_IDs(is_connected_cluster_small) = [];
            N_connected_clusters = numel(Connected_cluster_IDs);
        end
        
        % Mark the neighbours found as visited.
        Neighbour_lin_inds = Neighbour_inds(:, 1)+(Neighbour_inds(:, 2)-1)*im_size(1);
        is_pixel_notvisited(Neighbour_lin_inds) = false;
        pixels_clusterID(Neighbour_lin_inds) = Neighbours_cluster_ID;
        
        % Remove neighbours' indices from the the regional pixel indices array.
        Region_pixels_global_inds = Region_pixels_global_inds(~is_region_pixel_in_neighborhood, :);
        
        % Iterate through each connected cluster to scoop regional
        % pixels.
        parent_cluster_position = reshape(Clusters_position(i-1, j, :), [1, 2]);
        parent_cluster_size = Clusters_size(i-1, j);
        for k = 1:N_connected_clusters
            Cluster_ID = Cluster_ID + 1;
            cluster_iter_ind = N_new_clusters + k;
            connected_cluster_pixel_inds = Neighbour_inds(Neighbours_cluster_ID == Connected_cluster_IDs(k), :);
            connected_cluster_pixel_positions = connected_cluster_pixel_inds - 0.5;
            
            % Calculate the center of mass of the connected cluster.
            if use_pixel_mass
                connected_cluster_pixel_lin_inds = connected_cluster_pixel_inds(:, 1)+(connected_cluster_pixel_inds(:, 2)-1)*im_size(1);
                connected_cluster_pixel_mass = options.Pixel_mass(connected_cluster_pixel_lin_inds);
                Total_mass = sum(connected_cluster_pixel_mass);
                cluster_com = sum(connected_cluster_pixel_positions.*repmat(connected_cluster_pixel_mass, 1, 2), 1)/Total_mass;
                if any(isnan(cluster_com))
                    error('The center of mass is undefined.');
                end
            else
                cluster_com = mean(connected_cluster_pixel_positions, 1);
            end
            
            % Calculate the size of the cluster from its bounding box.
            rowcol_minmax_temp = [min(connected_cluster_pixel_inds, [], 1);max(connected_cluster_pixel_inds, [], 1)];
            cluster_size = sqrt(sum((diff(rowcol_minmax_temp)+1).^2, 2));
            
            % Calculate the smallest ratio between the cluster size and its
            % parent cluster size.
            if parent_cluster_size <= cluster_size
                size_ratio = parent_cluster_size/cluster_size;
            else
                size_ratio = cluster_size/parent_cluster_size;
            end
            
            % Calculate the position of the cluster.
            cluster_position = parent_cluster_position + (0.5)^size_ratio*(cluster_com-parent_cluster_position);
            %cluster_position = parent_cluster_position + 0.25*(cluster_com-parent_cluster_position);
            
            % Calculate the scooping distance. This distance is the maximal
            % distance between a cluster point and its position.
            scooping_dist_squared = max(sum(bsxfun(@minus, connected_cluster_pixel_positions, cluster_position).^2, 2));
            
            % Find all unvisited pixels in the parent cluster region that are within
            % a scooping distance of the cluster position. Add these pixels
            % to the connected cluster pixels to form the cluster neighborhood pixels.
            Region_pixels_positions = Region_pixels_global_inds - 0.5;
            region_pixel_dist_squared = sum(bsxfun(@minus, Region_pixels_positions, cluster_position).^2, 2);
            is_region_pixel_within_cluster = region_pixel_dist_squared <= scooping_dist_squared;
            additional_pixels_inds = Region_pixels_global_inds(is_region_pixel_within_cluster, :);
            cluster_neighborhood_pixels_inds = [connected_cluster_pixel_inds; additional_pixels_inds];
            
            % Among the cluster neighborhood pixels, retain only pixels
            % which are connected to the cluster.
            cluster_neighborhood_pixels_inds_ID = find_connected_clusters(cluster_neighborhood_pixels_inds);
            cluster_pixels_inds = cluster_neighborhood_pixels_inds(cluster_neighborhood_pixels_inds_ID==cluster_neighborhood_pixels_inds_ID(1),:);
            
            % Remove region pixels so that they are not scooped by other
            % connected clusters.
            Region_pixels_global_inds = Region_pixels_global_inds(~is_region_pixel_within_cluster, :);
            
            % Record cluster's pixels inds, position and size.
            Clusters_position(i, cluster_iter_ind, :) = cluster_position;
            Clusters_size(i, cluster_iter_ind) = cluster_size;
            Clusters_pixel_inds{i, cluster_iter_ind} = cluster_pixels_inds;
            Clusters_ID(i, cluster_iter_ind) = Cluster_ID;
            Clusters_parentID(i, cluster_iter_ind) = Clusters_ID(i-1, j);
            
            % Mark cluster pixels as visited.
            for l = 1:size(cluster_pixels_inds)
                is_pixel_notvisited(cluster_pixels_inds(l, 1), cluster_pixels_inds(l, 2)) = false;
                pixels_clusterID(cluster_pixels_inds(l, 1), cluster_pixels_inds(l, 2)) = j*k;
            end
        end
        
        % Increment the total number of clusters for the next iteration by
        % the number of connected clusters found for this parent cluster (j).
        N_new_clusters = N_new_clusters + N_connected_clusters;
        
    end
    
    if i >= 24 && options.PlotIterations && mod(i-1, 1) == 0
        % Plot a colored image of the initial binary image including the
        % clusters's found in this iteration and the previous iteration.
        image_to_plot = double(pixels_clusterID_previous > 0) + image_bin + (pixels_clusterID+1).*double(pixels_clusterID > 0 & pixels_clusterID_previous == 0);
        Iterations_f = figure;
        Iterations_f.Position(3:4) = 800;
        
        imagesc(image_to_plot);
        colorbar;
        axis equal;
        
        % Plot nodes positions.
        N_clusters_plot = N_clusters(i);
        clusters_pos = reshape(Clusters_position(i, 1:N_clusters_plot, :), [], 2);
        hold on
        plot(clusters_pos(:, 2), clusters_pos(:, 1), 'wx')
        
        close(Iterations_f);
    end
    
    % Move to next iteration.
    N_clusters(i) = N_new_clusters;
    i = i+1;
end

if options.Waitbar
    close(wb);
end

% Remove clusters that were not defined.
%is_cluster_defined = ~isnan(Clusters_ID);
is_cluster_defined = Clusters_ID > 0;
Clusters_position = reshape(Clusters_position, [], 2);
Clusters_position = Clusters_position(is_cluster_defined, :);
Clusters_ID = Clusters_ID(is_cluster_defined);
Clusters_parentID = Clusters_parentID(is_cluster_defined);
%% Calculate the local branch radius at each cluster position.
if options.Waitbar
    wb = waitbar(0,'Calculating the branch thickness at each cluster position.');
    Clusters_branch_thickness = find_branch_thickness(Clusters_position, image_bin);
    waitbar(1,wb); close(wb);
else
    Clusters_branch_thickness = find_branch_thickness(Clusters_position, image_bin);
end
%% Change coordinate system.
% Clusters_position is in in pixel coordinates. Define a new position in
% the  usual X,Y coordinates with the origin at the bottom left position.
Clusters_position_XY = Clusters_position(:, [2, 1]);
Clusters_position_XY(:, 2) = im_size(1)-(Clusters_position_XY(:, 2)-1);
%% Construct the Tree structure from the nodes structure.
% Define nodes properties from the clusters' position. Nodes are scanned to
% build the tree branches.
Nodes_Position = Clusters_position_XY;
Nodes_Diameter = Clusters_branch_thickness;
Nodes_ParentID = Clusters_parentID;
Nodes_ID = Clusters_ID;
N_nodes = size(Nodes_Position, 1);

% Define the position of the Soma as the first node position.
SomaPos = Nodes_Position(1,:);

% Redefine the Nodes_ID without gaps.
New_Nodes_ID = (1:N_nodes)';
Nodes_ParentID = rep(Nodes_ParentID, Nodes_ID, New_Nodes_ID);
if any(~ismember(Nodes_ParentID(2:end), New_Nodes_ID))
    error('Some parents ID were not replaced properly');
end
Nodes_ID = New_Nodes_ID;

% Determine the type of the node based on the Parent ID.
% Count the number of times a node is the parent of another node.
% Based on this count, the nodes will be defined as:
% 0:endpoint, 1:core node, >=2:branchpoint
ParIDscounts = histc(Nodes_ParentID, Nodes_ID);
Nodes_Type = ParIDscounts;
Nodes_Type(Nodes_Type > 1) = 2;

% Skeletonize the tree by visiting all nodes from endpoints to the root node (ParentID = 0).
Endpoints_NodeID = find(Nodes_Type == 0);
N_endpoints = numel(Endpoints_NodeID);
N_branchpoints = nnz(Nodes_Type == 2 & Nodes_ParentID >= 1);
Nodes_BranchID = nan(1, N_nodes);

N_Branches = N_endpoints+N_branchpoints;
Tree = struct();
Tree(N_Branches).PointsPos = [];
Tree(N_Branches).PointsDiameter = [];
Tree(N_Branches).Points_ClusterID = [];

BranchID = 0;
for n = 1:N_endpoints
    NodeID = Endpoints_NodeID(n);
    BranchID = BranchID+1;
    
    % Trace until an unconnected node or a node that has been
    % visited before is encountered.
    node_wasvisited = false;
    while NodeID >= 1 && ~node_wasvisited
        % Add node to the current branch.
        Tree(BranchID).Points_ClusterID = [NodeID; Tree(BranchID).Points_ClusterID];
        Tree(BranchID).PointsPos = [Nodes_Position(NodeID, :); Tree(BranchID).PointsPos];
        Tree(BranchID).PointsDiameter = [Nodes_Diameter(NodeID, :); Tree(BranchID).PointsDiameter];
        
        % If a branchpoint (NodeType == 2) is reached or if a root node (ParentID < 1) is
        % reached, terminate the current branch.
        NodeType = Nodes_Type(NodeID);
        ParentnodeID = Nodes_ParentID(NodeID);
        node_wasvisited = ~isnan(Nodes_BranchID(NodeID));
        if NodeType == 2 || ParentnodeID < 1
            % Assign the branch ID of the node.
            if ParentnodeID < 1
                Nodes_BranchID(NodeID) = 0;
            elseif ~node_wasvisited
                Nodes_BranchID(NodeID) = BranchID+1;
            end
            
            % Assign the parent ID of the current branch.
            Tree(BranchID).ParentID = Nodes_BranchID(NodeID);
            
            % Start a new branch if the parent node is connected (>=1) and the node was not visited.
            if ParentnodeID >= 1 && ~node_wasvisited
                BranchID = BranchID+1;
                Tree(BranchID).Points_ClusterID = [NodeID; Tree(BranchID).Points_ClusterID];
                Tree(BranchID).PointsPos = Nodes_Position(NodeID, :);
                Tree(BranchID).PointsDiameter = Nodes_Diameter(NodeID, :);
            end
        else
            % If the node is not a branchpoint, record its branch ID.
            Nodes_BranchID(NodeID) = BranchID;
        end
        
        % Move to the parent node.
        NodeID = ParentnodeID;
    end
end
%% Define the length and dimensionalfull length length_dim.
for n = 1:N_Branches
    Tree(n).Length = size(Tree(n).PointsPos, 1)-1;
    if Tree(n).Length > 0
        Tree(n).Length_dim = sum(sqrt(sum(diff(Tree(n).PointsPos).^2, 2)));
    else
        Tree(n).Length_dim = 0;
    end
end
%% Assign the parent, children and sibling ID for each branch.
Tree(N_Branches).ChildrenID = [];

% Define children ID.
for n = 1:N_Branches
    ParID = Tree(n).ParentID;
    if ParID ~= 0
        Tree(ParID).ChildrenID = [n, Tree(ParID).ChildrenID];
    end
end

% Define sibling ID
for n = 1:N_Branches
    ChildrenIDs = Tree(n).ChildrenID;
    if ~isempty(ChildrenIDs)
        for j = 1:numel(ChildrenIDs)
            Tree(ChildrenIDs(j)).SiblingID = ChildrenIDs([1:j-1, j+1:end]);
        end
    end
end

% Define sibling IDs for initial branches.
InitialBranches = find([Tree.ParentID] == 0);
if numel(InitialBranches) >= 2
    for i = 1:numel(InitialBranches)
        Tree(InitialBranches(i)).SiblingID = setdiff(InitialBranches, InitialBranches(i));
    end
else
    Tree(InitialBranches).SiblingID = 0;
end

% Order the tree.
Tree = order_tree(Tree);
%% Prune branches.
% Remove branches that have a path length smaller than the minimum. 
if isa(options.L_min, 'double') && options.L_min > 0
    ShortbranchesID = find([Tree.Length_dim] < options.L_min & cellfun(@isempty, {Tree.ChildrenID}));
    while ~isempty(ShortbranchesID)
        Tree = delete_branches(Tree, ShortbranchesID);
        N_Branches = numel(Tree);
        for n = 1:N_Branches
            % Recalculate the length.
            Tree(n).Length = size(Tree(n).PointsPos, 1)-1;
            Tree(n).Length_dim = sum(sqrt(sum(diff(Tree(n).PointsPos).^2, 2)));
        end
        
        ShortbranchesID = find([Tree.Length_dim] < options.L_min & cellfun(@isempty, {Tree.ChildrenID}));
    end
elseif isa(options.L_min, 'char') && strcmp(options.L_min, 'local_diameter')
    % Define the branch length threshold based on the local branch
    % diameter.
    while 1
        N_Branches = numel(Tree);
        for n = 1:N_Branches
            ParID = Tree(n).ParentID;
            if ParID > 0
                % Calculate the average branch diameter along the points of the
                % branch and use that as a length threshold.
                Tree(n).Length_thres = mean(Tree(ParID).PointsDiameter(:));
            else
                Tree(n).Length_thres = Tree(n).PointsDiameter(1);
            end
            Tree(n).Length_dim = sum(sqrt(sum(diff(Tree(n).PointsPos).^2, 2)));
        end
        
        ShortbranchesID = find([Tree.Length_dim] < [Tree.Length_thres] & cellfun(@isempty, {Tree.ChildrenID}));
        
        % Break of the loop if there are no more branches to remove.
        if isempty(ShortbranchesID)
            break;
        end
        
        % Prune.
        Tree = delete_branches(Tree, ShortbranchesID);
    end
end
%% Remove unecessary Tree fields and order the remaining fields.
Tree = rmfield(Tree, 'Points_ClusterID');

if isfield(Tree,'Length_thres')
    Tree = rmfield(Tree, 'Length_thres');
end
Tree = orderfields(Tree,{'PointsPos','PointsDiameter','Length','Length_dim','ParentID','SiblingID','ChildrenID','Depth'});
%% Center Tree at the Soma.
PositionOffset = options.SomaPos - SomaPos;
for i = 1:numel(Tree)
    Tree(i).PointsPos = Tree(i).PointsPos + PositionOffset;
end
%% Resample the branches' nodes at a constant rate.
Tree = resample_branches(Tree, options.Branch_nodes_distance);
%% Plot the final tree structure overlaid on top of the image.
if options.Plots && exist('plottree',2)
    %plottree(Tree, 'Lengthscale', 1, 'BackgroundImage', {image_input, 1, Soma_PixelInd});
end
end