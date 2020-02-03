function [Tree,DisconnectedBranchIDs] = disconnect_branches(Tree,DisconnectedBranchIDs)
% Function that disconnects input branches and connects the
% respective sibling branches to their parents. For branches with
% children branches, children are connected to the parent branch.

% The disconnected branch IDs are given as output, but not deleted from the tree
% structure.

% Input
% Tree = structure containing the tree branches.

% Output
% DisconnectedBranchIDs = IDs of the branches that will be disconnected.
%% Determine fieldnames for children and sibling branches.
Tree_fieldnames = fieldnames(Tree);

% Determine the field name of the children IDs.
if isfield(Tree,'DaughtersID')
    ChildrenID_name = 'DaughtersID';
elseif isfield(Tree,'ChildrenID')
    ChildrenID_name = 'ChildrenID';
else
    error('The children IDs fieldname cannot be determined');
end

% Determine the field name of the sibling IDs.
if isfield(Tree,'SisterID')
    SiblingID_name = 'SisterID';
elseif isfield(Tree,'SiblingID')
    SiblingID_name = 'SiblingID';
else
    error('The children IDs fieldname cannot be determined');
end

% Determine the field name of the branches' nodes.
if isfield(Tree,'PointsInd')
    Points_name = 'PointsInd';
elseif isfield(Tree,'PointsPos')
    Points_name = 'PointsPos';
else
    error('The children IDs fieldname cannot be determined');
end

% Find other array fields that have the same number of rows as the Points field.
% These field arrays will also be rearranged.
Field_size = cellfun(@(x) size(x),struct2cell(Tree(1)),'Uni',0);
Field_NDims = cellfun(@numel,Field_size);
Array_fieldnames = Tree_fieldnames(cellfun(@(x) x(1),Field_size) == size(Tree(1).(Points_name),1) & Field_NDims <= 2);

% Remove children ID fieldname from the array fieldnames.
Array_fieldnames = Array_fieldnames(~ismember(Array_fieldnames,ChildrenID_name));
N_array_fieldnames = numel(Array_fieldnames);
%% Disconnect branches
% For each branchID that will be deleted, connect its sibling to its parent.
N_DisconnectedBranches = numel(DisconnectedBranchIDs);
AdditionalDisconnectedBranchIDs = nan(1,N_DisconnectedBranches);
for i=1:N_DisconnectedBranches
    BranchID = DisconnectedBranchIDs(i);
    
    % Get the parent ID.
    ParentID = Tree(BranchID).ParentID;
    ChildrenIDs = Tree(BranchID).(ChildrenID_name);
    SiblingIDs = Tree(BranchID).(SiblingID_name);
    
    % Verify that the branch is connected. If the branch is disconnected,
    % move on to the next branch. To verify if a branch is connected,
    % verify that ID references are consistent. For example, verify that
    % the Parent ID of the branch is among the children IDs of the parent.
    if ParentID > 0
        is_ParentID_consistent = any(Tree(ParentID).(ChildrenID_name) == BranchID);
    else
        is_ParentID_consistent = ~isnan(ParentID);
    end
    
    if any(SiblingIDs > 0)
        is_SiblingID_consistent = all(cellfun(@(x) any(x==BranchID),{Tree(SiblingIDs(SiblingIDs>0)).(SiblingID_name)}));
    else
        is_SiblingID_consistent = all(~isnan(SiblingIDs));
    end
    
    if all(~isnan(ChildrenIDs))
        is_ChildrenID_consistent = all([Tree(ChildrenIDs).ParentID] == BranchID);
    else
        is_ChildrenID_consistent = 0;
    end
    
    is_branch_connected = is_ParentID_consistent & is_SiblingID_consistent & is_ChildrenID_consistent;
    
    if ~is_branch_connected
        continue;
    end

    % Update the branch references.
    % Update the children's branches (if any).
    for j = ChildrenIDs
        % Update the parent ID.
        Tree(j).ParentID = ParentID;
        
        % Update the points position.
        %tree(j).PointsPos(1,:) = tree(BranchID).PointsPos(1,:);
        for k=1:N_array_fieldnames
            Tree(j).(Array_fieldnames{k})(1,:) = Tree(BranchID).(Array_fieldnames{k})(1,:);
        end
    end
    
    % Remove the branchID from the parent's ChildrenID.
    parent_new_ChildrenIDs = [Tree(BranchID).(SiblingID_name) ChildrenIDs];
    parent_new_ChildrenIDs = parent_new_ChildrenIDs(parent_new_ChildrenIDs>0);
    if ParentID > 0
        %parent_new_ChildrenIDs = [tree(ParentID).(ChildrenID_name) ChildrenIDs];
        %parent_new_ChildrenIDs = parent_new_ChildrenIDs(parent_new_ChildrenIDs~=BranchID);
        Tree(ParentID).(ChildrenID_name) = parent_new_ChildrenIDs;
    end
    
    % Reassign the SiblingID of the new parent's children.
    for j=1:numel(parent_new_ChildrenIDs)
        ID = parent_new_ChildrenIDs(j);
        Tree(ID).(SiblingID_name) = parent_new_ChildrenIDs(parent_new_ChildrenIDs~=ID);
    end
    
    % After branch deletion, if the parent has a single child and it is not
    % being deleted, join the child and parent's branches together by
    % making the parent branch inherit the child's branch nodes.
    N_children = numel(parent_new_ChildrenIDs);
    if N_children == 1 && parent_new_ChildrenIDs > 0 && ~ismember(parent_new_ChildrenIDs,DisconnectedBranchIDs)
        ChildID = parent_new_ChildrenIDs;
        
        % Add the length and points of the child branch to the parent branch.
        ParentLength = Tree(ParentID).Length;
        SiblingLength = double(Tree(ChildID).Length);
        Tree(ParentID).Length = ParentLength+SiblingLength;
        
        % Modify the field arrays.
        for k=1:N_array_fieldnames
            Tree(ParentID).(Array_fieldnames{k})(ParentLength+1:ParentLength+SiblingLength+1, :) = Tree(ChildID).(Array_fieldnames{k})(1:SiblingLength+1, :);
        end
        Tree(ParentID).(ChildrenID_name) = Tree(ChildID).(ChildrenID_name);
        
        % Update the parent ID of the child's children.
        Child_ChildrenIDs = Tree(ChildID).(ChildrenID_name);
        for j = Child_ChildrenIDs
            Tree(j).ParentID = ParentID;
        end
        
        % Delete the child branch since its content has been added to
        % the parent.
        AdditionalDisconnectedBranchIDs(i) = ChildID;
    end
    
    % Remove ID references in the disconnected branch.
    Tree(BranchID).ParentID = nan;
    Tree(BranchID).(ChildrenID_name) = nan;
    Tree(BranchID).(SiblingID_name) = nan;
end

%% Output
AdditionalDisconnectedBranchIDs = AdditionalDisconnectedBranchIDs(~isnan(AdditionalDisconnectedBranchIDs));
DisconnectedBranchIDs = [DisconnectedBranchIDs AdditionalDisconnectedBranchIDs];