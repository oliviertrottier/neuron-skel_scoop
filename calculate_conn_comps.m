function [N_connected_comp,levels] = calculate_conn_comps(im,varargin)
% Function to find the number of connected components above a set of
% intensity threshold levels.

% Input
% im = image (2D array)

% Output
% N_connected_comp = number of connected components.

% Calling options
%
% calculate_conn_comps(im)
% Calculate the number of connected components at 10 equally-spaced threshold levels
% between the minium and maximum pixel intensity.
%
% calculate_conn_comps(im,N)
% Calculate the number of connected components at N equally-spaced threshold levels
% between the minium and maximum pixel intensity.
%
% calculate_conn_comps(im,levels)
% Calculate the number of connected components at each input threshold 'levels'.
%%
im = double(im);
if nargin > 1
    if isdouble(varargin{1}) && numel(varargin{1})==1
        N_levels = varargin{1};
        levels = linspace(min(im(:)),max(im(:)),N_levels+2);
    elseif isdouble(varargin{1})
        levels = [min(im(:)) varargin{1} max(im(:))];
        N_levels = numel(levels) - 2;
    else
        N_levels = 10;
        levels = linspace(min(im(:)),max(im(:)),N_levels+2);
    end
end

N_connected_comp = zeros(N_levels,1);
levels = levels(2:end-1);
component_min_size = 10;
for i=1:N_levels
    im_bin = im > levels(i);
    conn_comps = bwconncomp(im_bin);
    N_connected_comp(i) = nnz(cellfun(@(x) numel(x) >= component_min_size,conn_comps.PixelIdxList));
end
end