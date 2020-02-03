function [neighbour_inds,is_candidate_in_neighborhood] = find_cluster_neigborhood(cluster_pixels_inds,candidate_pixels_inds,neigborhood_radius)
% Function to find which candidate pixels are in the neigborhood of cluster
% pixels.

% Input
% cluster_pixels_inds = cluster pixel 2D indices.
% candidate_pixels_inds = candidate pixel 2D indices.
% neigborhood_radius = radius of the neighborhood used for search.

% Output
% neighbour_inds = indices of the candidate pixels that are neighbours of
%                  the cluster pixels.
% is_candidate_in_neighborhood = logical array indicating if the input
%                               candidate pixel are in the neighborhood.
%%
% Calculate the smallest chebychev distance between each candidate pixel and a cluster pixel.
candidate_pixel_dist = pdist2(cluster_pixels_inds,candidate_pixels_inds,'chebychev','Smallest',1);
%candidate_pixel_dist = pdist2(cluster_pixels_inds,candidate_pixels_inds,'euclidean','Smallest',1);

% Determine which candidate pixels are in the neighborhood. If the
% chebyshev distance is used, the search mask is a squared. If the
% euclidean distance is used, the search mask is a circle.
is_candidate_in_neighborhood = candidate_pixel_dist(:) <= neigborhood_radius;
neighbour_inds = candidate_pixels_inds(is_candidate_in_neighborhood,:);
end