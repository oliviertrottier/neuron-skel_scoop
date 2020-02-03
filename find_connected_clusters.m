function [pixels_cluster_ID,N_clusters,Clusters_size] = find_connected_clusters(pixel_indices)
% Function to find the connected clusters of a set of pixel indices.

% Input
% pixel_indices = N x N_d array consisting of the indices of N pixels in
% N_d dimensions.

% Output
% pixels_cluster_ID = N x 1 array giving the cluster ID of each pixel.
% N_clusters = number of clusters found.
% Clusters_size = size of each cluster (nnz(pixels_cluster_ID==i) = Clusters_size(i)).

[N_pixels,N_dimensions] = size(pixel_indices);
N_neighbours_max = 3*N_dimensions - 1;

is_pixel_visited = false(N_pixels,1);
pixels_cluster_ID = nan(N_pixels,1);
pixels_connected_IDs = zeros(N_pixels,N_neighbours_max);
Clusters_size = nan(10,1);

% Store the ID of pixels used to probe other connected pixels. The ID
% correspond to the row index of pixel_indices. 
pixel_probes_ID = 1;
N_clusters = 1;
Cluster_size = 1;
N_unvisited_pixels = N_pixels - 1;
is_pixel_visited(pixel_probes_ID) = true;
pixels_cluster_ID(pixel_probes_ID) = N_clusters;

while N_unvisited_pixels > 0
    % Get the index of the current pixel.
    N_probes = numel(pixel_probes_ID);
    
    % Find all unvisited pixels that are neigbours of each probe.
    if N_probes > 0
        pixel_probe_inds = pixel_indices(pixel_probes_ID,:);
        for i=1:N_probes
            pixel_prob_ind = pixel_probe_inds(i,:);
            Unvisited_pixel_IDs = find(~is_pixel_visited);
            Unvisited_pixel_inds = pixel_indices(Unvisited_pixel_IDs,:);
            
            % Find unvisited pixels that are connected to pixel_ind.
            Unvisited_pixel_dist = bsxfun(@minus,Unvisited_pixel_inds,pixel_prob_ind);
            is_pixel_connected = all(abs(Unvisited_pixel_dist)<=1,2);
            connected_pixels_ID_temp = Unvisited_pixel_IDs(is_pixel_connected);
            N_connected_pixels = numel(connected_pixels_ID_temp);
            
            if N_connected_pixels > 0
                % Record the pixels IDs connected to pixel_IDs(i).
                pixels_connected_IDs(pixel_probes_ID(i),1:N_connected_pixels) = connected_pixels_ID_temp;
                
                % Record the cluster ID of the connected pixels.
                pixels_cluster_ID(connected_pixels_ID_temp) = N_clusters;
                
                % Indicate the connected pixels as visited.
                is_pixel_visited(connected_pixels_ID_temp) = true;
                
                % Decrement the number of unvisited pixels.
                N_unvisited_pixels = N_unvisited_pixels - N_connected_pixels;
                
                % Increment the size of the cluster.
                Cluster_size = Cluster_size + N_connected_pixels;
            end
        end
        
        % Define the newly-connected pixels as the probes for the
        % next iteration.
        pixel_probes_ID = nonzeros(pixels_connected_IDs(pixel_probes_ID,:));
    else
        % If no more probes exist for the current cluster, start a new
        % cluster.
        % Record the cluster size.
        Clusters_size(N_clusters) = Cluster_size;
        
        % Find another unvisited pixel and start a new cluster.
        pixel_probes_ID = find(~is_pixel_visited,1);
        N_unvisited_pixels = N_unvisited_pixels - 1;
        N_clusters = N_clusters + 1;
        Cluster_size = 1;
        is_pixel_visited(pixel_probes_ID) = true;
        pixels_cluster_ID(pixel_probes_ID) = N_clusters;
    end
end
Clusters_size(N_clusters) = Cluster_size;
Clusters_size = Clusters_size(1:N_clusters);

% Send error if some pixel have no cluster IDs.
if any(isnan(pixels_cluster_ID))
    error('Not all connected components were found.');
end
end