function Tree_ordered = order_tree(Tree)
% Function to order a Tree branches in ascending depth and branching angle.
% Root branches (depth = 1), are defined as branches whose Parent ID=0.

% Input
% Tree = Input tree structure that.

% Output
% Tree_ordered = Order tree structure.
%% Determine field names of sibling and children,
% Determine the field name of the children IDs.
if isfield(Tree, 'DaughtersID')
    ChildrenID_name = 'DaughtersID';
elseif isfield(Tree, 'ChildrenID')
    ChildrenID_name = 'ChildrenID';
else
    error('The children IDs fieldname cannot be determined');
end

% Determine the field name of the sibling IDs.
if isfield(Tree, 'SisterID')
    SiblingID_name = 'SisterID';
elseif isfield(Tree, 'SiblingID')
    SiblingID_name = 'SiblingID';
else
    error('The children IDs fieldname cannot be determined');
end

if ~isfield(Tree,'PointsPos')
    error('PointsPos is needed to order the Trees');
end
%% Define the depth of each branch.
% Find the number of branches at the soma and initialize their depth to 1.
ParentIDs = reshape([Tree.ParentID],[],1);
BranchesID = find(ParentIDs==0);

depth = 1;
while ~isempty(BranchesID)
    % Assign depth.
    [Tree(BranchesID).Depth] = deal(depth);
    
    % Move to the next depth level.
    BranchesID = [Tree(BranchesID).(ChildrenID_name)];
    depth = depth + 1;
end
Depths = [Tree(:).Depth]';
%% Calculate the branching angle of each branch. 
% This will be used to order the branches.
NBranches = numel(Tree);
BranchingAngles = nan(NBranches, 1);
if isfield(Tree,'PointsPos')
    for i = 1:NBranches
        if size(Tree(i).PointsPos,1) > 1
            vec = diff(Tree(i).PointsPos(1:2, :));
            BranchingAngles(i) = mod(atan2(vec(2), vec(1)),2*pi);
        end
    end
end
%% Calculate the orientation of each branch's branchpoint.
Branchpoint_thetas = nan(NBranches, 1);
Soma_position = Tree(1).PointsPos(1, :);
for i=1:NBranches
    vec = Tree(i).PointsPos(1, :) - Soma_position;
    Branchpoint_thetas(i) = atan2(vec(2), vec(1));
end
%% Reorder the Tree in ascending depth, branchpoint orientation and branching angle.
% First, find the map that maps the old branch IDs to the new IDs.
% i = OldOrder(NewOrder(i))
% NewOrder(i) gives the old branch Id of the branch whose new ID is i.
[~,NewOrder] = sortrows([Depths Branchpoint_thetas BranchingAngles]);
% OldOrder(i) gives the new branch ID of the branch whose old ID is i.
[~,OldOrder] = sort(NewOrder);

Tree_ordered = reshape(Tree(NewOrder),size(Tree));
for i=1:NBranches
    ParID = Tree_ordered(i).ParentID;
    if ParID > 0
        Tree_ordered(i).ParentID = OldOrder(ParID);
    end
    
    SiblingIDs = Tree_ordered(i).(SiblingID_name);
    if SiblingIDs > 0
        Tree_ordered(i).(SiblingID_name) = OldOrder(SiblingIDs)';
    end
    
    ChildrenIDs = Tree_ordered(i).(ChildrenID_name);
    Tree_ordered(i).(ChildrenID_name) = OldOrder(ChildrenIDs)';
end
end