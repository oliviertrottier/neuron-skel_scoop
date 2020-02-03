function [im_filled,im_edges] = im_gradient_fill(im)
% Function that finds edges in an image and fills its interior by following the intensity gradient.

% First, the edges of the shape are found using the Canny method. Then,
% tracers are initialized at each edge pixel and moved along the image
% intensity gradient with a small amount of rotational diffusion.

% Input
% im = 2D array representing the image.

% Output
% im_filled = 2D logical array representing the interior of the shapes
%             found in the input image.
% im_edges = 2D logical array representing the edges found in the input
%            image.
%%
% Calculate the image gradient.
[grad_mag,grad_dir] = imgradient(im);
grad_mag_normalized = grad_mag/max(grad_mag(:)); % This is what is thresholded in the edge() function.
angular_noise_sigma = 30;

% Find pixels that belong to edges in the input image.
mag_thres = graythresh(grad_mag_normalized);
canny_thresholds = [0.4 1]*mag_thres;
%canny_thresholds = multithresh(grad_mag_normalized,2);
[im_edges,~] = edge(im,'Canny',canny_thresholds);
im_size = size(im);

% Remove edge clusters that are smaller in size than a given threshold.
min_edge_cluster_size = 10;
edge_conn_comp = bwconncomp(im_edges,8);
cluster_sizes = cellfun(@numel,edge_conn_comp.PixelIdxList);
small_clusters_logind = cluster_sizes <= min_edge_cluster_size;
small_clusters_linind = cell2mat(edge_conn_comp.PixelIdxList(small_clusters_logind)');
im_edges(small_clusters_linind) = false;

% Test plots
% [~,biggest_cluster_ind] = max(cluster_sizes);
% im_label = bwlabel(im_edges,8);
% figure;imshowpair(im,im_label==biggest_cluster_ind);
% figure;histogram(im(im_edges));
% figure;imshowpair(im,im_edges_filled);

% figure;histogram(im);
% 
% bin_edges = 0:double(max(im(im_edges)));
% counts = histc(im(im_edges),bin_edges);
% otsu_thres = otsuthresh(counts);
% bin_edges(round(otsu_thres*numel(bin_edges)))
% 
% bin_edges = 0:double(max(im(:)));
% counts = histc(im(:),bin_edges);
% otsu_thres = otsuthresh(counts);
% bin_edges(round(otsu_thres*numel(bin_edges)))

% %% Annihilation phase. 
% % Move edge pixels for a certain number of
% % displacements and remove them if they encounter an initial edge pixel.
% N_displacements_annihilation = 10;
% is_edge_pixel_removed = false(nnz(im_edges),1);
% [edge_pix_row,edge_pix_col] = find(im_edges);
% edge_pix_lin = edge_pix_row + (edge_pix_col-1)*im_size(1);
% for i=1:N_displacements_annihilation
%     % Calculate the displacement of each edge pixel and use it to calculate their new position.
%     %pix_displacement_angle = grad_dir_rounded(edge_pix_lin);
%     pix_displacement_angle = round((grad_dir(edge_pix_lin))/45)*45;
%     
%     pix_ind_displacement = sign([-sind(pix_displacement_angle) cosd(pix_displacement_angle)]);
%     new_pix_ind = [edge_pix_row,edge_pix_col] + pix_ind_displacement;
%     
%     % Make sure the indices are not out of bound.
%     new_pix_ind(:,1) = max(min(new_pix_ind(:,1),im_size(1)),1);
%     new_pix_ind(:,2) = max(min(new_pix_ind(:,2),im_size(2)),1);
%     
%     % Define the new edge pixel linear index.
%     new_pix_lin_ind = new_pix_ind(:,1) + (new_pix_ind(:,2)-1)*im_size(1);
%     
%     % Remove the pixel if its new pixel position is an original edge pixel.
%     is_edge_pixel_removed = is_edge_pixel_removed | im_edges(new_pix_lin_ind);
%     
%     % Prepare for next iteration.
%     edge_pix_row = new_pix_ind(:,1);
%     edge_pix_col = new_pix_ind(:,2);
%     edge_pix_lin = new_pix_lin_ind;
% end
% 
% % Remove the edge pixels that collided with original edge pixels.
% %im_edges_temp = im_edges;
% [edge_pix_row,edge_pix_col] = find(im_edges);
% edge_pix_lin = edge_pix_row + (edge_pix_col-1)*im_size(1);
% im_edges(edge_pix_lin(is_edge_pixel_removed)) = false;
%% Fill the interior of the image by walking tracers along the intensity gradients.
% Move the edge pixels along the image gradient to fill the interior of the image.
% Each pixel is moved a given number of times (N_displacements). The process is repeated for a
% certain number of iterations (N_iterations).
N_iterations = 5;
N_displacements = 500;
im_filled = im_edges;
for j=1:N_iterations
    im_edges_temp = im_edges;
    for i=1:N_displacements
        [edge_pix_row,edge_pix_col] = find(im_edges_temp);
        edge_pix_lin = edge_pix_row + (edge_pix_col-1)*im_size(1);
        
        % Calculate the displacement of each edge pixel and use it to calculate their new position.
        %pix_displacement_angle = grad_dir_rounded(edge_pix_lin);
        pix_displacement_angle = round((grad_dir(edge_pix_lin) + angular_noise_sigma*randn(size(edge_pix_lin)))/45)*45;
        
        pix_ind_displacement = sign([-sind(pix_displacement_angle) cosd(pix_displacement_angle)]);
        new_pix_ind = [edge_pix_row,edge_pix_col] + pix_ind_displacement;
        
        % Make sure the indices are not out of bound.
        new_pix_ind(:,1) = max(min(new_pix_ind(:,1),im_size(1)),1);
        new_pix_ind(:,2) = max(min(new_pix_ind(:,2),im_size(2)),1);
        
        % Define a new edge image with the displaced pixels.
        new_pix_lin_ind = new_pix_ind(:,1) + (new_pix_ind(:,2)-1)*im_size(1);
        new_im_edges = false(im_size);
        new_im_edges(new_pix_lin_ind) = true;
        
        % Save the position of the displaced edge pixels in the filled
        % image.
        im_filled(new_pix_lin_ind) = true;
        
        % Prepare for next iteration.
        im_edges_temp = new_im_edges;
    end
end

% Fill in holes that are smaller than 4 pixels.
holes_max_size = 4;
im_holes = ~im_filled;
holes_conn_comp = bwconncomp(im_holes,4);
holes_sizes = cellfun(@numel,holes_conn_comp.PixelIdxList);
small_holes_logind = holes_sizes <= holes_max_size;
small_holes_linind = cell2mat(holes_conn_comp.PixelIdxList(small_holes_logind)');
im_filled(small_holes_linind) = true;

end