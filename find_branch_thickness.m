function [thickness, theta_perp, is_cross] = find_branch_thickness(pixel_indices, image_intensity, varargin)
% Function to find the branch radius and perpendicular angle at given
% positions in an image.

% Input
% pixel_indices = pixel indices where the thickness will be computed.
%                 Each row corresponds to a different position.
%                 If the input array is an n x 2 array, subscript indices are assumed.
%                 If the input array is an n x 1 array, linear indices are assumed.
%                 Non-integer positions are rounded.
% image_intensity = pixel intensity of the image.

% Output
% thickness = branch thickness at each input position.
% theta_perp = angle of the direction perpendicular to the branch at each input position.
% is_cross = logical array indicating if the input pixel position is on an
%            edge.
%% Parse input
image_size = size(image_intensity);
switch size(pixel_indices,2)
    case 1
        % Check that the linear indices are not out of bound.
        assert(all(pixel_indices>=1) & all(pixel_indices<=prod(image_size)),'Pixel indices are out of bounds');
        
        % Convert the pixel indices to subscripted indices.
        pixel_indices = [mod(pixel_indices - 1,image_size(1))+1 ceil(pixel_indices/image_size(1))];
    case 2
        % Check that the subscript indices are within the range of the
        % image size.
        assert(all(pixel_indices(:,1)>=1) &...
            all(pixel_indices(:,1)<=image_size(1)) &...
            all(pixel_indices(:,2)>=1) &...
            all(pixel_indices(:,2)<=image_size(2)),'Pixel indices are out of bounds');
    otherwise
        error('Pixel indices must be a N x 1 or N x 2 array');
end
%% Parse optional parameters
p = inputParser;
% Statistic to compute with the thickness of all line segments at each given point.
addParameter(p, 'Statistic', 'min', @(x) ismember(x, {'min', 'max', 'mean'}));
addParameter(p, 'Plots', false); % Plot intermediate plots detailing the method.
addParameter(p, 'Detect_cross',false); % Detect branch crosses.
parse(p, varargin{:});
options = p.Results;
%% Define the pixel coordinates along each line segments.
statistic_func = str2func(options.Statistic);

r_max = 100;
dr = .01*r_max;
rs = 0:dr:r_max;
N_rs = numel(rs);
r_mid_ind = ceil(numel(rs)/2);
r_mid = rs(r_mid_ind);
d_theta = 0.01*pi;
thetas = 0:d_theta:(pi-d_theta);
thetas = [thetas, thetas+pi];
N_thetas = numel(thetas);
d_theta = diff(thetas(1:2));
[theta_grid, r_grid] = ndgrid(thetas, rs);
r_grid_pos = r_grid(:, 2:end);

% Angular average range for edge pixels.
% For pixels on an edge, average the thickness of lines that are within an
% angle 'edge_theta_search_range' from the mean angle of lines directed
% inward.
edge_theta_search_range = pi/8;
edge_theta_search_range_ind = round(-edge_theta_search_range/d_theta):round(edge_theta_search_range/d_theta);

% Define the pixel coordinates of the line at different angles. In order
% for the angle to correspond to the usual cartesian angle, we need to
% transform theta to -theta. This way, at theta=pi/2, for example, the line row coordinates
% (y = row, x = col) will be negative for positive r(corresponding to upward offset in the image array);
lines_pixel_coordinates = cat(3, r_grid.*sin(-theta_grid), r_grid.*cos(-theta_grid));
%% Upsample the image by interpolating the original image.
pixel_row_ind = (1:image_size(1))-0.5;
pixel_col_ind = (1:image_size(2))-0.5;
image_gridded_interpolant = griddedInterpolant({pixel_row_ind, pixel_col_ind}, double(image_intensity));

% Upsample with the interpolant.
Upsampling_factor = 1;
pixel_row_ind_fine = (pixel_row_ind(1)-1+1/Upsampling_factor):1/Upsampling_factor:(pixel_row_ind(end));
pixel_col_ind_fine = (pixel_col_ind(1)-1+1/Upsampling_factor):1/Upsampling_factor:pixel_col_ind(end);
[pixels_fine_row_ind, pixels_fine_col_ind] = ndgrid(pixel_row_ind_fine, pixel_col_ind_fine);
image_intensity_fine = image_gridded_interpolant(pixels_fine_row_ind, pixels_fine_col_ind);
image_fine_size = size(image_intensity_fine);

% Find the new lines' indices in terms of the upsampled image.
lines_pixel_coordinates_fine = round(lines_pixel_coordinates*Upsampling_factor);
pixel_positions_fine = round(pixel_indices*Upsampling_factor);

% Determine the intensity threshold that will be used to determine if a
% pixel is part of a branch.
intensity_threshold = mean(image_intensity_fine(:)) + 0.5*std(image_intensity_fine(:));
max_thickness = 2*(r_max-dr);

N_positions = size(pixel_positions_fine, 1);
thickness = nan(N_positions, 1);
theta_perp = nan(N_positions, 1);
if options.Detect_cross
    is_cross = false(N_positions,1);
else
    is_cross = [];
end
%% Scan through all pixel positions and find its thickness.
pixel_positions_fine = shiftdim(pixel_positions_fine, -1);
for i = 1:N_positions
    % For each position, calculate the interpolated intensity along half-lines
    % that originate from each pixel position at different angle.
    lines_coordinates_temp = bsxfun(@plus, lines_pixel_coordinates_fine, pixel_positions_fine(1, i, :));
    
    % Remove out of range indices.
    lines_row_ind = lines_coordinates_temp(:, :, 1);
    lines_row_ind(lines_row_ind < 1) = 1;
    lines_row_ind(lines_row_ind > image_fine_size(1)) = image_fine_size(1);
    lines_coordinates_temp(:, :, 1) = lines_row_ind;
    
    lines_col_ind = lines_coordinates_temp(:, :, 2);
    lines_col_ind(lines_col_ind < 1) = 1;
    lines_col_ind(lines_col_ind > image_fine_size(2)) = image_fine_size(2);
    lines_coordinates_temp(:, :, 2) = lines_col_ind;
    
    % Change to linear indices.
    lines_coordinates_temp_linind = lines_coordinates_temp(:, :, 1) + image_fine_size(1)*(lines_coordinates_temp(:, :, 2)-1);
    lines_intensities = image_intensity_fine(lines_coordinates_temp_linind);
    
    % Find the thickness of the branch projected along the 
    % different scanning lines (which have different angles).
    
    % To find the projected thickness along a certain angle,
    % find the radial positions where the intensity falls below the threshold.
    lines_intensities_thres = lines_intensities > intensity_threshold;
    %intensity_change = [ones(N_thetas, 1), diff(lines_intensities_thres, 1, 2)];
    intensity_change = diff(lines_intensities_thres, 1, 2);
    intensity_drop_ind = intensity_change == -1;
    intensity_drop_ind(:,end) = 1;
    %intensity_drop_rs = rs(end)*ones(N_thetas, N_rs-1);
    intensity_drop_rs = nan(N_thetas, N_rs-1);
    intensity_drop_rs(intensity_drop_ind) = r_grid_pos(intensity_drop_ind);
    [~,intensity_drop_smallest_r_ind] = min(intensity_drop_rs, [], 2);
    lines_length = rs(intensity_drop_smallest_r_ind+1)'-dr;
    
    % Calculate the projected thicknesses by adding the length of a line with
    % its symmetrically opposite line.
    Projected_thicknesses = lines_length(1:N_thetas/2) + lines_length(N_thetas/2+1:end);
    lines_ind = 1:N_thetas/2;
    
    % Determine if the query position is on an edge. The pixel is on an
    % edge if any one of the half line have a zero length.
    is_line_zero_length = lines_length == 0;
    is_on_edge = any(is_line_zero_length);
    
    % If the pixel is on an edge, find the line segments that are
    % one-sided, i.e., segments whose thresholded intensities is only in the
    % negative or positive radial positions. Only these segments should be
    % considered in the thickness calculation.
    if is_on_edge
        % Find the average angle of the zero-length lines.
        zero_length_line_theta = thetas(is_line_zero_length)';
        zero_length_line_vec_avg = mean([cos(zero_length_line_theta), sin(zero_length_line_theta)]);
        zero_length_line_angle_avg = mod(atan2(zero_length_line_vec_avg(2), zero_length_line_vec_avg(1)), pi);
        
        % Find the indices of the lines that within the angular search
        % range of the mean angle.
        zero_length_line_angle_avg_ind = round(zero_length_line_angle_avg/d_theta);
        line_segments_candidate_ind = mod(zero_length_line_angle_avg_ind + edge_theta_search_range_ind - 1, N_thetas/2)+1;
        
        % Do not remove lines that are within the average angle.
        is_line_removed = true(N_thetas/2, 1);
        is_line_removed(line_segments_candidate_ind) = false;
        
        % Remove line segments.
        Projected_thicknesses(is_line_removed) = [];
        lines_ind(is_line_removed) = [];
        
    end
    
    % Compute the thickness and the perpendicular angle.
    thickness(i) = statistic_func(Projected_thicknesses);
    [thickness_min, min_thickness_ind] = min(Projected_thicknesses);
    has_line_thickness_min = Projected_thicknesses == thickness_min;
    
    % If there is more than one line that has the minimal thickness,
    % choose the line whose angle is closest to the average line segment
    % angle. 
    if nnz(has_line_thickness_min) > 1
        first_min_cluster_end_ind = find(diff(has_line_thickness_min)==-1,1);
        if isempty(first_min_cluster_end_ind)
            min_thickness_ind = round(mean(lines_ind(has_line_thickness_min)));
        else
            min_thickness_ind = round(mean(lines_ind(has_line_thickness_min(1:first_min_cluster_end_ind))));
        end
    end
    theta_perp(i) = thetas(min_thickness_ind);
    
    % Determine if the point is a potential cross by finding the number of local
    % length maxima along the angular direction. If the position is a
    % cross, there will be 4 local maxima in the line length.
    if options.Detect_cross
        scan_length = 1.5*thickness(i);
        if nnz(lines_length > 1.5*thickness(i)) > 4
            lines_length_temp = [lines_length(min_thickness_ind+1:end); lines_length(1:min_thickness_ind)];
            
            %[~,Peaks_loc] = findpeaks(lines_length_temp,'MinPeakHeight',1.5*thickness(i),'MinPeakProminence',10*dr,'MinPeakDistance',round(pi/8/d_theta));
            isabove_thresh = lines_length_temp >= scan_length;
            thres_change = diff(isabove_thresh);
            start_peak_ind = find(thres_change==1) + 1;
            end_peak_ind = find(thres_change==-1);
            
            % Check if there are exactly 4 peaks in the length of line
            % segments.
            if numel(start_peak_ind) == 4 && numel(end_peak_ind) == 4
                Peaks_loc = mean([start_peak_ind end_peak_ind],2);
                Peaks_width = end_peak_ind - start_peak_ind;
                
                % If there are 4 peaks, check that the peaks are diametrically opposite.
                % In other words, check that the distance between a peak
                % and its second-next neighbor is ~180 deg.
                Peak_loc_diff = [Peaks_loc(3) - Peaks_loc(1) Peaks_loc(4) - Peaks_loc(2)];
                peak_dist_min = round((pi - pi/8)/d_theta);
                peak_dist_max = round((pi + pi/8)/d_theta);
                is_loc_diff_180 = Peak_loc_diff >= peak_dist_min & Peak_loc_diff <= peak_dist_max;
                are_peak_diametrically_opposite = all(is_loc_diff_180);
                
                if are_peak_diametrically_opposite
                    % Next, check that the width of diametrically opposite
                    % peaks are similar. We expect this to be true at a cross
                    % section.
                    opposite_peak_width_diff_max = 2;
                    Peaks_width_diff = abs([Peaks_width(3) - Peaks_width(1) Peaks_width(4) - Peaks_width(2)]);
                    is_cross(i) = all(Peaks_width_diff < opposite_peak_width_diff_max);
                end
            end
        end
    end
    
    % Send warning if calculate thickness is above the maximum.
    if thickness(i) >= max_thickness
        warning('The calculated thickness (%.2f) is above the maximum thickness (%.2f). The thickness may be inaccurate', thickness(i), max_thickness);
    end
    %% Plot the intensity along each line.
    if options.Plots
        % Transform angles to degrees.
        thetas_d = thetas/pi*180;
        
        % Plot the image intensities near the query position.
        neigh_size = r_max;
        center_ind = pixel_positions_fine(1, i, :);
        neigh_row_ind = max(-neigh_size+center_ind(1), 1):min(neigh_size+center_ind(1), image_size(1));
        neigh_col_ind = max(-neigh_size+center_ind(2), 1):min(neigh_size+center_ind(2), image_size(2));
        image_neigh = double(image_intensity(neigh_row_ind, neigh_col_ind));
        
        query_pixel_ind = [find(neigh_row_ind == center_ind(1)), find(neigh_col_ind == center_ind(2))];
        image_neigh(query_pixel_ind(1), query_pixel_ind(2)) = 2;
        figure;
        imagesc(image_neigh)
        
        % Plot the line intensities at each r,theta.
        figure;
        imagesc([rs(1), rs(end)], [thetas_d(1), thetas_d(end)], lines_intensities);
        hold on
        plot([rs(1), rs(end)], ones(1, 2)*thetas_d(min_thickness_ind), 'w--', 'DisplayName', sprintf('Thickness = %.2f, \\theta = %.2f ^\\circ', thickness(i), thetas_d(min_thickness_ind)));
        xlabel('Radius');
        ylabel('\theta (^\circ)');
        
        
        % Plot the edge angle, if such a case happen.
        if is_on_edge
            % Calculate the average angle of the removed line segments.
            removed_segments_angle = thetas(is_line_removed)';
            removed_segments_vec_avg = mean([cos(removed_segments_angle), sin(removed_segments_angle)]);
            removed_segments_angle_avg = atan2(removed_segments_vec_avg(2), removed_segments_vec_avg(1));
            
            removed_segments_angle_avg_d = removed_segments_angle_avg/pi*180;
            plot([rs(1), rs(end)], ones(1, 2)*removed_segments_angle_avg_d, 'r--', 'DisplayName', sprintf('Removed angles average, \\theta = %.2f ^\\circ', removed_segments_angle_avg_d));
            
            % Plot the omitted lines segments.
            line_segments_removed_ind = find(is_line_removed);
            patch_break_ind = find([0;diff(line_segments_removed_ind) ~= 1]);
            
            patch_alpha = 0.3;
            if isempty(patch_break_ind)
                % Plot the single patch
                corners_r = [rs(1), rs(end), rs(end), rs(1)];
                corners_theta = [ones(1, 2)*thetas_d(line_segments_removed_ind(1)), ones(1, 2)*thetas_d(line_segments_removed_ind(end))];
                patch(corners_r, corners_theta, 'red', 'FaceAlpha', patch_alpha, 'DisplayName', 'OmittedSegments')
            else
                % Plot two patches.
                patch1_corners_r = [rs(1), rs(end), rs(end), rs(1)];
                patch1_corners_theta = [ones(1, 2)*thetas_d(line_segments_removed_ind(1)), ones(1, 2)*thetas_d(line_segments_removed_ind(patch_break_ind-1))];
                patch(patch1_corners_r, patch1_corners_theta, 'red', 'FaceAlpha', patch_alpha, 'DisplayName', 'OmittedSegments')
                
                patch2_corners_r = [rs(1), rs(end), rs(end), rs(1)];
                patch2_corners_theta = [ones(1, 2)*thetas_d(line_segments_removed_ind(patch_break_ind)), ones(1, 2)*thetas_d(line_segments_removed_ind(end))];
                patch(patch2_corners_r, patch2_corners_theta, 'red', 'FaceAlpha', patch_alpha, 'DisplayName', 'OmittedSegments')
            end
        end
        
        % Plot legend
        legend;
    end
end
end