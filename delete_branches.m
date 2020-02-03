function Tree = delete_branches(Tree, DeletedBranchIDs, varargin)
%% Parse optional parameters
p = inputParser;
% Disconnect branches before deletion. 
% Disconnection removes references from and to a branch, but doesn't remove it from the Tree structure.
addParameter(p, 'Disconnect', true); 
parse(p, varargin{:});
options = p.Results;
options_isdefault = cell2struct(num2cell(ismember(p.Parameters, p.UsingDefaults)), p.Parameters, 2);
%%
% Function that deletes branches with ID=branchIDs and connects the
% respective sister branches to their parents.
tree_fieldnames = fieldnames(Tree);

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

% First, disconnect branches before deleting them.
if options.Disconnect
    [Tree, BranchIDs_Removed] = disconnect_branches(Tree, DeletedBranchIDs);
else
    BranchIDs_Removed = DeletedBranchIDs;
end
%% Assign branch ID
% Find an ID name to use for storing the Branch ID.
IsIDNameExistent = true;
% Record the rng state to reverse it after finding a fieldname.
Rng=rng();
while IsIDNameExistent
    IDname = ['ID', num2str(randi(1000))];
    IsIDNameExistent = ismember(IDname, tree_fieldnames);
end
rng(Rng);

% Assign ID to be used for further reference.
IDs = num2cell(1:numel(Tree));
[Tree.(IDname)] = deal(IDs{:});
IDs = cell2mat(IDs);
BranchIDs_Kept = setdiff(IDs,BranchIDs_Removed);
%% Perform sanity checks on deleted and non-deleted branches.
% Check if deleted branches are referenced by non-deleted branches.
Tree_new = Tree(BranchIDs_Kept);
ParentIDs_Kept = [Tree_new.ParentID];
if any(ismember(ParentIDs_Kept, BranchIDs_Removed))
    error('Some removed branches are referenced as Parent IDs.')
end
SiblingIDs_Kept = [Tree_new.(SiblingID_name)];
if any(ismember(SiblingIDs_Kept, BranchIDs_Removed))
    error('Some removed branches are referenced as Sibling IDs.')
end
ChildrenIDs_Kept = [Tree_new.(ChildrenID_name)];
if any(ismember(ChildrenIDs_Kept, BranchIDs_Removed))
    error('Some removed branches are referenced as Children IDs.')
end

% Check that the Sibling ID of deleted branchs is non-empty.
SiblingIDs = {Tree_new.(SiblingID_name)};
has_empty_siblings = cellfun(@isempty, SiblingIDs);
if any(has_empty_siblings)
    has_empty_siblings_ID = BranchIDs_Kept(has_empty_siblings);
    error('Some branches have no Sibling ID.')
end

% Check that the Parent ID of deleted branchs is non-empty.
ParentIDs = {Tree_new.ParentID};
has_empty_parent = cellfun(@isempty, ParentIDs);
if any(has_empty_parent)
    has_empty_parent_ID = BranchIDs_Kept(has_empty_parent);
    error('Some branches have no Parent ID.')
end

%% Delete the branches from the tree structure.
% Remove the branches.
% isdeleted = ismember(IDs, BranchIDs_Removed);
% Tree_new = Tree(~isdeleted);
Tree = Tree_new;
%% Fix ID references before reordering the branches.
OldIDs = [Tree.(IDname)];
N_Branches = numel(Tree);
NewIDs = 1:N_Branches;

ParentIDs = rep([Tree.ParentID], OldIDs, NewIDs);
SiblingsID = rep({Tree.(SiblingID_name)}, OldIDs, NewIDs);
ChildrenIDs = rep({Tree.(ChildrenID_name)}, OldIDs, NewIDs);
for i = 1:N_Branches
    Tree(i).(IDname) = NewIDs(i);
    Tree(i).ParentID = ParentIDs(i);
    Tree(i).(SiblingID_name) = SiblingsID{i};
    Tree(i).(ChildrenID_name) = ChildrenIDs{i};
end
%% Reorder the tree structure.
Tree = order_tree(Tree);

% Check if ID references go out of bounds.
allIDs = cell2mat(cellfun(@double,{[Tree.ParentID], [Tree.(SiblingID_name)], [Tree.(ChildrenID_name)]},'Uni',0));
if any(allIDs > N_Branches)
    error('Some ID replacements are incorrect.')
end

% Remove ID fieldname.
Tree = rmfield(Tree, IDname);
end