function Interpolated_line_positions = interp_equi(Points_position,Segment_length)
% Interpolate N-D points to form a line. The segments of the line have a
% length of Segment_length.
% Points are assumed to be ordered.

% Input
% Points_position = N x N_d array representing the N_d-D positions of N points.
% Segment_length = length of each segment in the final line.

% Output
% Interpolated_line_positions = N x N_d array representing the points that outline the
%               interpolated line.
%%
[N_points,N_dim] = size(Points_position);
assert(N_points >= 2,'At least 2 points are necessary.')

Interpolated_line_positions = zeros(N_points,N_dim);
N_interp_nodes = 1;
Segment_length_squared = Segment_length^2;
distance_threshold = Segment_length_squared/2;

% Find the starting position.
% Find all points within a segment length of the first points and average
% their position to define the starting position.
squared_dists = sum(bsxfun(@minus,Points_position,Points_position(1,:)).^2,2);

segment_start_pos = mean(Points_position(squared_dists < Segment_length_squared,:),1);
Interpolated_line_positions(1,:) = segment_start_pos;
segment_start_point_ind = 1;
Points_ind = 2;

while Points_ind <= N_points
    dist_squared = sum((Points_position(Points_ind,:) - segment_start_pos).^2,2);
    if dist_squared >= distance_threshold
        % Fit all points between the start of the segment to the current
        % point position.
        Points_to_fit = Points_position(segment_start_point_ind:Points_ind,:);
        
        % Find the direction of the segment.
        % If there are many points contained in the current segment, find the best fit line by performing PCA on the points.
        % If there are only two points, direct the segment towards the
        % second point.
        if size(Points_to_fit,1) > 2
            [Pca_vecs,Points_pca_score] = pca(Points_to_fit);
            
            % Determine the orientation by checking the pca score of the last
            % fitted points. If the score is negative, it implies that the
            % segment is oriented away from the last point. In this case, flip
            % the orientation of the segment 180 deg.
            Segment_orientation = 2*double(Points_pca_score(end,1) > 0) - 1;
            Segment_unit_vec = Segment_orientation*Pca_vecs(:,1)';
        else
            Segment_unit_vec = Points_to_fit(2,:) - segment_start_pos;
            Segment_unit_vec = Segment_unit_vec/sqrt(sum(Segment_unit_vec.^2,2));
        end
        segment_end_pos = segment_start_pos + Segment_length*Segment_unit_vec;
        
        % When interpolating the last point, stop the interpolation if the
        % segment's end position is moving away from the last point.
        if Points_ind == N_points
            segment_start_dist = sum(segment_start_pos - Points_position(N_points,:).^2,2);
            segment_end_dist = sum(segment_end_pos - Points_position(N_points,:).^2,2);
            if segment_start_dist < segment_end_dist
                break;
            end
        end
        
        % Save the newly interpolated point.
        N_interp_nodes = N_interp_nodes + 1;
        Interpolated_line_positions(N_interp_nodes,:) = segment_end_pos;
        
        % Start the next segment.
        segment_start_pos = segment_end_pos;
        segment_start_point_ind = Points_ind - 1;
    else
        % Go to the next point.
        Points_ind = Points_ind + 1;
    end
end

% Truncate unecessary nodes in the interpolated line.
Interpolated_line_positions = Interpolated_line_positions(1:N_interp_nodes,:);

% Plot the interpolated line if no output is requested.
if nargout==0
    figure;
    plot(Points_position(:,1),Points_position(:,2),'x','DisplayName','Input')
    axis equal;
    hold on;
    plot(Interpolated_line_positions(:,1),Interpolated_line_positions(:,2),'o','DisplayName','Interpolation')
    legend;
end
end