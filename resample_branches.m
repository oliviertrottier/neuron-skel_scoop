function Tree = resample_branches(Tree,stepsize)
% Function to resample the branches of a tree structure so that branches
% are made of nodes that are a distance "stepsize" apart.

% Input
% Tree = structure containing the tree branches.
% stepsize = stepsize between each node of each branch.
%%
N_branches = numel(Tree);

% Use the entry in Tree to find all fieldnames that the same number of rows
% as PointsPos. These fieldnames will also be resampled.
Tree_fieldnames = fieldnames(Tree);
Fieldnames_size = cellfun(@(x) size(x),struct2cell(Tree(1)),'Uni',0);
Fieldnames_NDims = cellfun(@numel,Fieldnames_size);

Fieldnames_to_resample = Tree_fieldnames(cellfun(@(x) x(1),Fieldnames_size) == size(Tree(1).PointsPos,1) & Fieldnames_NDims <= 2);
Fieldnames_to_resample = setdiff(Fieldnames_to_resample,'PointsPos');
N_Fieldnames_to_resample = numel(Fieldnames_to_resample);
Field_sizes = zeros(1,N_Fieldnames_to_resample+1);
for i=1:N_branches
    % Change the branchpoint position to the parent endpoint. This is
    % necessary because the parent branch endpoint may have been moved in a
    % previous iteration of the interpolation.
    ParID = Tree(i).ParentID;
    if ParID > 0
        Tree(i).PointsPos(1,:) = Tree(ParID).PointsPos(end,:);
        
        % Change other fields.
        for k=1:N_Fieldnames_to_resample
            Tree(i).(Fieldnames_to_resample{k})(1,:) = Tree(ParID).(Fieldnames_to_resample{k})(end,:);
        end
    end
    
    % Resample the points positions.
    PointsPos = Tree(i).PointsPos;
    if size(PointsPos,1) < 2
        error('Cannot resample with only one point.')
    end
    
    % Calculate the step sizes for the current branch nodes and keep only
    % the non-zero stepsizes.
    stepsizes = sqrt(sum(diff(PointsPos).^2,2));
    is_stepsize_nonzero = stepsizes > 1e-5;
    sample_log_ind = [true;is_stepsize_nonzero];
    stepsizes = [0; stepsizes(is_stepsize_nonzero)];
    PointsPos = PointsPos(sample_log_ind,:);
    
    % Interpolate if there is more than 1 point on the branch after removing zero stepsizes.
    if size(PointsPos,1) > 1
        Length_cum = cumsum(stepsizes);
        NewLength_cum = (0:stepsize:(Length_cum(end)+stepsize))';
        NewPointsPos = interp1(Length_cum, PointsPos, NewLength_cum,'linear','extrap');
        
        % Test that the points were interpolated correctly.
        New_stepsizes = sqrt(sum(diff(NewPointsPos).^2,2));
        are_steps_even = all(abs(New_stepsizes - stepsize) < 1e4*eps);
        if ~are_steps_even
            % Interpolate with a user-define interpolation that ensures
            % equal step sizes.
            NewPointsPos = interp_equi(PointsPos,stepsize);
            
%             % Interpolate with a spline at half the step size and
%             % reinterpolate linearly.
%             Dense_Length_cum = (0:stepsize/2:(Length_cum(end) + stepsize/2))';
%             Dense_PointsPos = interp1(Length_cum, PointsPos, Dense_Length_cum,'spline','extrap');
%             
%             Dense_stepsizes = sqrt(sum(diff(Dense_PointsPos).^2,2));
%             Dense_Length_cum = [0;cumsum(Dense_stepsizes)];
%             
%             NewPointsPos = interp1(Dense_Length_cum, Dense_PointsPos, NewLength_cum,'linear','extrap');
        end
        
        % Find which endpoint to use between the last two new points (if
        % there is more than 2 points).
        if size(NewPointsPos,1) > 2
            PointsPos_Cand = NewPointsPos(end-1:end,:);
            PointsPos_Cand = PointsPos_Cand - PointsPos(end,:);
            PointsPos_Cand_dist2 = sum(PointsPos_Cand.^2,2);
            [~,bestcand_ind] = min(PointsPos_Cand_dist2);
        else
            bestcand_ind = 2;
        end
        NewPointsPos = NewPointsPos(1:end-2+bestcand_ind,:);
        
        % Interpolate the other fields.
        for k=1:numel(Fieldnames_to_resample)
            if are_steps_even
                Field_interpolated = interp1(Length_cum, Tree(i).(Fieldnames_to_resample{k})(sample_log_ind,:), NewLength_cum,'linear','extrap');
                Tree(i).(Fieldnames_to_resample{k}) = Field_interpolated(1:end-2+bestcand_ind,:);
            else
                Field_val = Tree(i).(Fieldnames_to_resample{k})(sample_log_ind,:);
                Interpolant = scatteredInterpolant(PointsPos(:,1),PointsPos(:,2),Field_val);
                Tree(i).(Fieldnames_to_resample{k}) = Interpolant(NewPointsPos(:,1),NewPointsPos(:,2));
            end
        end
    else
        NewPointsPos = PointsPos;
        for k=1:numel(Fieldnames_to_resample)
            Tree(i).(Fieldnames_to_resample{k}) = Tree(i).(Fieldnames_to_resample{k})(1,:);
        end
    end
    
    % Assing the new points position.
    Tree(i).PointsPos = NewPointsPos;
    
    % Update the length of the branch.
    Tree(i).Length = size(NewPointsPos,1)-1;
    
    % Check that all resampled fields have the same number of rows.
    Field_sizes(1) = size(Tree(i).PointsPos,1);
    for k=1:numel(Fieldnames_to_resample)
        Field_sizes(k+1) = size(Tree(i).(Fieldnames_to_resample{k}),1);
    end
    assert(numel(unique(Field_sizes))==1,'The Points field are not of equal size');
end

% % Plot the interpolated line.
% if 0
%     figure;
%     plot(PointsPos(:,1),PointsPos(:,2),'x')
%     hold on;plot(NewPointsPos(:,1),NewPointsPos(:,2),'o')
%     NewPointsPos2 = interp_equi(PointsPos,stepsize);
%     plot(NewPointsPos2(:,1),NewPointsPos2(:,2),'s');
%     
%     figure;
%     dist = cumsum([0; sqrt(sum(diff(PointsPos,1).^2,2))]);
%     newdist = cumsum([0; sqrt(sum(diff(NewPointsPos,1).^2,2))]);
%     plot(dist,'x')
%     hold on;plot(newdist,'o')
% end
end